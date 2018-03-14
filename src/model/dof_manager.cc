/**
 * @file   dof_manager.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Aug 18 2015
 * @date last modification: Wed Feb 21 2018
 *
 * @brief  Implementation of the common parts of the DOFManagers
 *
 * @section LICENSE
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "dof_manager.hh"
#include "communicator.hh"
#include "element_group.hh"
#include "mesh.hh"
#include "mesh_utils.hh"
#include "node_group.hh"
#include "non_linear_solver.hh"
#include "sparse_matrix.hh"
#include "time_step_solver.hh"
/* -------------------------------------------------------------------------- */
#include <memory>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
DOFManager::DOFManager(const ID & id, const MemoryID & memory_id)
    : Memory(id, memory_id),
      communicator(Communicator::getStaticCommunicator()) {}

/* -------------------------------------------------------------------------- */
DOFManager::DOFManager(Mesh & mesh, const ID & id, const MemoryID & memory_id)
    : Memory(id, memory_id), mesh(&mesh), local_system_size(0),
      pure_local_system_size(0), system_size(0),
      communicator(mesh.getCommunicator()) {
  this->mesh->registerEventHandler(*this, _ehp_dof_manager);
}

/* -------------------------------------------------------------------------- */
DOFManager::~DOFManager() = default;

/* -------------------------------------------------------------------------- */
// void DOFManager::getEquationsNumbers(const ID &, Array<UInt> &) {
//   AKANTU_TO_IMPLEMENT();
// }

/* -------------------------------------------------------------------------- */
std::vector<ID> DOFManager::getDOFIDs() const {
  std::vector<ID> keys;
  for (const auto & dof_data : this->dofs)
    keys.push_back(dof_data.first);

  return keys;
}

/* -------------------------------------------------------------------------- */
void DOFManager::assembleElementalArrayLocalArray(
    const Array<Real> & elementary_vect, Array<Real> & array_assembeled,
    const ElementType & type, const GhostType & ghost_type, Real scale_factor,
    const Array<UInt> & filter_elements) {
  AKANTU_DEBUG_IN();

  UInt nb_element;
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_degree_of_freedom =
      elementary_vect.getNbComponent() / nb_nodes_per_element;

  UInt * filter_it = nullptr;
  if (filter_elements != empty_filter) {
    nb_element = filter_elements.size();
    filter_it = filter_elements.storage();
  } else {
    nb_element = this->mesh->getNbElement(type, ghost_type);
  }

  AKANTU_DEBUG_ASSERT(elementary_vect.size() == nb_element,
                      "The vector elementary_vect("
                          << elementary_vect.getID()
                          << ") has not the good size.");

  const Array<UInt> & connectivity =
      this->mesh->getConnectivity(type, ghost_type);

  Array<Real>::const_matrix_iterator elem_it =
      elementary_vect.begin(nb_degree_of_freedom, nb_nodes_per_element);

  for (UInt el = 0; el < nb_element; ++el, ++elem_it) {
    UInt element = el;
    if (filter_it != nullptr) {
      // conn_it = conn_begin + *filter_it;
      element = *filter_it;
    }

    // const Vector<UInt> & conn = *conn_it;
    const Matrix<Real> & elemental_val = *elem_it;
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      UInt offset_node = connectivity(element, n) * nb_degree_of_freedom;
      Vector<Real> assemble(array_assembeled.storage() + offset_node,
                            nb_degree_of_freedom);
      Vector<Real> elem_val = elemental_val(n);
      assemble.aXplusY(elem_val, scale_factor);
    }

    if (filter_it != nullptr)
      ++filter_it;
    //    else
    //      ++conn_it;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DOFManager::assembleElementalArrayToResidual(
    const ID & dof_id, const Array<Real> & elementary_vect,
    const ElementType & type, const GhostType & ghost_type, Real scale_factor,
    const Array<UInt> & filter_elements) {
  AKANTU_DEBUG_IN();

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_degree_of_freedom =
      elementary_vect.getNbComponent() / nb_nodes_per_element;
  Array<Real> array_localy_assembeled(this->mesh->getNbNodes(),
                                      nb_degree_of_freedom);

  array_localy_assembeled.clear();

  this->assembleElementalArrayLocalArray(
      elementary_vect, array_localy_assembeled, type, ghost_type, scale_factor,
      filter_elements);

  this->assembleToResidual(dof_id, array_localy_assembeled, 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DOFManager::assembleElementalArrayToLumpedMatrix(
    const ID & dof_id, const Array<Real> & elementary_vect,
    const ID & lumped_mtx, const ElementType & type,
    const GhostType & ghost_type, Real scale_factor,
    const Array<UInt> & filter_elements) {
  AKANTU_DEBUG_IN();

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_degree_of_freedom =
      elementary_vect.getNbComponent() / nb_nodes_per_element;
  Array<Real> array_localy_assembeled(this->mesh->getNbNodes(),
                                      nb_degree_of_freedom);

  array_localy_assembeled.clear();

  this->assembleElementalArrayLocalArray(
      elementary_vect, array_localy_assembeled, type, ghost_type, scale_factor,
      filter_elements);

  this->assembleToLumpedMatrix(dof_id, array_localy_assembeled, lumped_mtx, 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DOFManager::assembleMatMulDOFsToResidual(const ID & A_id,
                                              Real scale_factor) {
  for (auto & pair : this->dofs) {
    const auto & dof_id = pair.first;
    auto & dof_data = *pair.second;

    this->assembleMatMulVectToResidual(dof_id, A_id, *dof_data.dof,
                                       scale_factor);
  }
}

/* -------------------------------------------------------------------------- */
DOFManager::DOFData::DOFData(const ID & dof_id)
    : support_type(_dst_generic), group_support("__mesh__"), dof(nullptr),
      blocked_dofs(nullptr), increment(nullptr), previous(nullptr),
      solution(0, 1, dof_id + ":solution"),
      local_equation_number(0, 1, dof_id + ":local_equation_number") {}

/* -------------------------------------------------------------------------- */
DOFManager::DOFData::~DOFData() = default;

/* -------------------------------------------------------------------------- */
DOFManager::DOFData & DOFManager::getNewDOFData(const ID & dof_id) {
  auto it = this->dofs.find(dof_id);
  if (it != this->dofs.end()) {
    AKANTU_EXCEPTION("This dof array has already been registered");
  }

  std::unique_ptr<DOFData> dofs_storage = std::make_unique<DOFData>(dof_id);
  this->dofs[dof_id] = std::move(dofs_storage);
  return *dofs_storage;
}

/* -------------------------------------------------------------------------- */
void DOFManager::registerDOFsInternal(const ID & dof_id,
                                      Array<Real> & dofs_array) {
  DOFData & dofs_storage = this->getDOFData(dof_id);
  dofs_storage.dof = &dofs_array;

  UInt nb_local_dofs = 0;
  UInt nb_pure_local = 0;

  const DOFSupportType & support_type = dofs_storage.support_type;

  switch (support_type) {
  case _dst_nodal: {
    UInt nb_nodes = 0;
    const ID & group = dofs_storage.group_support;

    NodeGroup * node_group = nullptr;
    if (group == "__mesh__") {
      nb_nodes = this->mesh->getNbNodes();
    } else {
      node_group = &this->mesh->getElementGroup(group).getNodeGroup();
      nb_nodes = node_group->size();
    }

    nb_local_dofs = nb_nodes;
    AKANTU_DEBUG_ASSERT(
        dofs_array.size() == nb_local_dofs,
        "The array of dof is too shot to be associated to nodes.");

    for (UInt n = 0; n < nb_nodes; ++n) {
      UInt node = n;
      if (node_group)
        node = node_group->getNodes()(n);

      nb_pure_local += this->mesh->isLocalOrMasterNode(node) ? 1 : 0;
    }

    nb_pure_local *= dofs_array.getNbComponent();
    nb_local_dofs *= dofs_array.getNbComponent();
    break;
  }
  case _dst_generic: {
    nb_local_dofs = nb_pure_local =
        dofs_array.size() * dofs_array.getNbComponent();
    break;
  }
  default: { AKANTU_EXCEPTION("This type of dofs is not handled yet."); }
  }

  this->pure_local_system_size += nb_pure_local;
  this->local_system_size += nb_local_dofs;

  communicator.allReduce(nb_pure_local, SynchronizerOperation::_sum);

  this->system_size += nb_pure_local;
}

/* -------------------------------------------------------------------------- */
void DOFManager::registerDOFs(const ID & dof_id, Array<Real> & dofs_array,
                              const DOFSupportType & support_type) {
  DOFData & dofs_storage = this->getNewDOFData(dof_id);
  dofs_storage.support_type = support_type;

  this->registerDOFsInternal(dof_id, dofs_array);
}

/* -------------------------------------------------------------------------- */
void DOFManager::registerDOFs(const ID & dof_id, Array<Real> & dofs_array,
                              const ID & support_group) {
  DOFData & dofs_storage = this->getNewDOFData(dof_id);
  dofs_storage.support_type = _dst_nodal;
  dofs_storage.group_support = support_group;

  this->registerDOFsInternal(dof_id, dofs_array);
}

/* -------------------------------------------------------------------------- */
void DOFManager::registerDOFsPrevious(const ID & dof_id, Array<Real> & array) {
  DOFData & dof = this->getDOFData(dof_id);

  if (dof.previous != nullptr) {
    AKANTU_EXCEPTION("The previous dofs array for "
                     << dof_id << " has already been registered");
  }

  dof.previous = &array;
}

/* -------------------------------------------------------------------------- */
void DOFManager::registerDOFsIncrement(const ID & dof_id, Array<Real> & array) {
  DOFData & dof = this->getDOFData(dof_id);

  if (dof.increment != nullptr) {
    AKANTU_EXCEPTION("The dofs increment array for "
                     << dof_id << " has already been registered");
  }

  dof.increment = &array;
}

/* -------------------------------------------------------------------------- */
void DOFManager::registerDOFsDerivative(const ID & dof_id, UInt order,
                                        Array<Real> & dofs_derivative) {
  DOFData & dof = this->getDOFData(dof_id);
  std::vector<Array<Real> *> & derivatives = dof.dof_derivatives;

  if (derivatives.size() < order) {
    derivatives.resize(order, nullptr);
  } else {
    if (derivatives[order - 1] != nullptr) {
      AKANTU_EXCEPTION("The dof derivatives of order "
                       << order << " already been registered for this dof ("
                       << dof_id << ")");
    }
  }

  derivatives[order - 1] = &dofs_derivative;
}

/* -------------------------------------------------------------------------- */
void DOFManager::registerBlockedDOFs(const ID & dof_id,
                                     Array<bool> & blocked_dofs) {
  DOFData & dof = this->getDOFData(dof_id);

  if (dof.blocked_dofs != nullptr) {
    AKANTU_EXCEPTION("The blocked dofs array for "
                     << dof_id << " has already been registered");
  }

  dof.blocked_dofs = &blocked_dofs;
}

/* -------------------------------------------------------------------------- */
void DOFManager::splitSolutionPerDOFs() {
  auto it = this->dofs.begin();
  auto end = this->dofs.end();

  for (; it != end; ++it) {
    DOFData & dof_data = *it->second;
    dof_data.solution.resize(dof_data.dof->size() *
                             dof_data.dof->getNbComponent());
    this->getSolutionPerDOFs(it->first, dof_data.solution);
  }
}

/* -------------------------------------------------------------------------- */
SparseMatrix &
DOFManager::registerSparseMatrix(const ID & matrix_id,
                                 std::unique_ptr<SparseMatrix> & matrix) {
  SparseMatricesMap::const_iterator it = this->matrices.find(matrix_id);
  if (it != this->matrices.end()) {
    AKANTU_EXCEPTION("The matrix " << matrix_id << " already exists in "
                                   << this->id);
  }

  SparseMatrix & ret = *matrix;
  this->matrices[matrix_id] = std::move(matrix);
  return ret;
}

/* -------------------------------------------------------------------------- */
/// Get an instance of a new SparseMatrix
Array<Real> & DOFManager::getNewLumpedMatrix(const ID & id) {
  ID matrix_id = this->id + ":lumped_mtx:" + id;
  LumpedMatricesMap::const_iterator it = this->lumped_matrices.find(matrix_id);
  if (it != this->lumped_matrices.end()) {
    AKANTU_EXCEPTION("The lumped matrix " << matrix_id << " already exists in "
                                          << this->id);
  }

  auto mtx =
      std::make_unique<Array<Real>>(this->local_system_size, 1, matrix_id);
  this->lumped_matrices[matrix_id] = std::move(mtx);
  return *this->lumped_matrices[matrix_id];
}

/* -------------------------------------------------------------------------- */
NonLinearSolver & DOFManager::registerNonLinearSolver(
    const ID & non_linear_solver_id,
    std::unique_ptr<NonLinearSolver> & non_linear_solver) {
  NonLinearSolversMap::const_iterator it =
      this->non_linear_solvers.find(non_linear_solver_id);
  if (it != this->non_linear_solvers.end()) {
    AKANTU_EXCEPTION("The non linear solver " << non_linear_solver_id
                                              << " already exists in "
                                              << this->id);
  }

  NonLinearSolver & ret = *non_linear_solver;
  this->non_linear_solvers[non_linear_solver_id] = std::move(non_linear_solver);

  return ret;
}

/* -------------------------------------------------------------------------- */
TimeStepSolver & DOFManager::registerTimeStepSolver(
    const ID & time_step_solver_id,
    std::unique_ptr<TimeStepSolver> & time_step_solver) {
  TimeStepSolversMap::const_iterator it =
      this->time_step_solvers.find(time_step_solver_id);
  if (it != this->time_step_solvers.end()) {
    AKANTU_EXCEPTION("The non linear solver " << time_step_solver_id
                                              << " already exists in "
                                              << this->id);
  }

  TimeStepSolver & ret = *time_step_solver;
  this->time_step_solvers[time_step_solver_id] = std::move(time_step_solver);
  return ret;
}

/* -------------------------------------------------------------------------- */
SparseMatrix & DOFManager::getMatrix(const ID & id) {
  ID matrix_id = this->id + ":mtx:" + id;
  SparseMatricesMap::const_iterator it = this->matrices.find(matrix_id);
  if (it == this->matrices.end()) {
    AKANTU_SILENT_EXCEPTION("The matrix " << matrix_id << " does not exists in "
                                          << this->id);
  }

  return *(it->second);
}

/* -------------------------------------------------------------------------- */
bool DOFManager::hasMatrix(const ID & id) const {
  ID mtx_id = this->id + ":mtx:" + id;
  auto it = this->matrices.find(mtx_id);
  return it != this->matrices.end();
}

/* -------------------------------------------------------------------------- */
Array<Real> & DOFManager::getLumpedMatrix(const ID & id) {
  ID matrix_id = this->id + ":lumped_mtx:" + id;
  LumpedMatricesMap::const_iterator it = this->lumped_matrices.find(matrix_id);
  if (it == this->lumped_matrices.end()) {
    AKANTU_SILENT_EXCEPTION("The lumped matrix "
                            << matrix_id << " does not exists in " << this->id);
  }

  return *(it->second);
}

/* -------------------------------------------------------------------------- */
const Array<Real> & DOFManager::getLumpedMatrix(const ID & id) const {
  ID matrix_id = this->id + ":lumped_mtx:" + id;
  auto it = this->lumped_matrices.find(matrix_id);
  if (it == this->lumped_matrices.end()) {
    AKANTU_SILENT_EXCEPTION("The lumped matrix "
                            << matrix_id << " does not exists in " << this->id);
  }

  return *(it->second);
}

/* -------------------------------------------------------------------------- */
bool DOFManager::hasLumpedMatrix(const ID & id) const {
  ID mtx_id = this->id + ":lumped_mtx:" + id;
  auto it = this->lumped_matrices.find(mtx_id);
  return it != this->lumped_matrices.end();
}

/* -------------------------------------------------------------------------- */
NonLinearSolver & DOFManager::getNonLinearSolver(const ID & id) {
  ID non_linear_solver_id = this->id + ":nls:" + id;
  NonLinearSolversMap::const_iterator it =
      this->non_linear_solvers.find(non_linear_solver_id);
  if (it == this->non_linear_solvers.end()) {
    AKANTU_EXCEPTION("The non linear solver " << non_linear_solver_id
                                              << " does not exists in "
                                              << this->id);
  }

  return *(it->second);
}

/* -------------------------------------------------------------------------- */
bool DOFManager::hasNonLinearSolver(const ID & id) const {
  ID solver_id = this->id + ":nls:" + id;
  auto it = this->non_linear_solvers.find(solver_id);
  return it != this->non_linear_solvers.end();
}

/* -------------------------------------------------------------------------- */
TimeStepSolver & DOFManager::getTimeStepSolver(const ID & id) {
  ID time_step_solver_id = this->id + ":tss:" + id;
  TimeStepSolversMap::const_iterator it =
      this->time_step_solvers.find(time_step_solver_id);
  if (it == this->time_step_solvers.end()) {
    AKANTU_EXCEPTION("The non linear solver " << time_step_solver_id
                                              << " does not exists in "
                                              << this->id);
  }

  return *(it->second);
}

/* -------------------------------------------------------------------------- */
bool DOFManager::hasTimeStepSolver(const ID & solver_id) const {
  ID time_step_solver_id = this->id + ":tss:" + solver_id;
  auto it = this->time_step_solvers.find(time_step_solver_id);
  return it != this->time_step_solvers.end();
}

/* -------------------------------------------------------------------------- */
void DOFManager::savePreviousDOFs(const ID & dofs_id) {
  this->getPreviousDOFs(dofs_id).copy(this->getDOFs(dofs_id));
}

/* -------------------------------------------------------------------------- */
/* Mesh Events                                                                */
/* -------------------------------------------------------------------------- */
std::pair<UInt, UInt>
DOFManager::updateNodalDOFs(const ID & dof_id, const Array<UInt> & nodes_list) {
  auto & dof_data = this->getDOFData(dof_id);
  UInt nb_new_local_dofs = 0;
  UInt nb_new_pure_local = 0;

  nb_new_local_dofs = nodes_list.size();
  for (const auto & node : nodes_list) {
    nb_new_pure_local += this->mesh->isLocalOrMasterNode(node) ? 1 : 0;
  }

  const auto & dof_array = *dof_data.dof;
  nb_new_pure_local *= dof_array.getNbComponent();
  nb_new_local_dofs *= dof_array.getNbComponent();

  this->pure_local_system_size += nb_new_pure_local;
  this->local_system_size += nb_new_local_dofs;

  UInt nb_new_global = nb_new_pure_local;
  communicator.allReduce(nb_new_global, SynchronizerOperation::_sum);

  this->system_size += nb_new_global;

  dof_data.solution.resize(dof_data.solution.size() + nb_new_local_dofs);

  return std::make_pair(nb_new_local_dofs, nb_new_pure_local);
}

/* -------------------------------------------------------------------------- */
void DOFManager::onNodesAdded(const Array<UInt> & nodes_list,
                              const NewNodesEvent &) {
  for (auto & pair : this->dofs) {
    const auto & dof_id = pair.first;
    auto & dof_data = this->getDOFData(dof_id);
    if (dof_data.support_type != _dst_nodal)
      continue;

    const auto & group = dof_data.group_support;

    if (group == "__mesh__") {
      this->updateNodalDOFs(dof_id, nodes_list);
    } else {
      const auto & node_group =
          this->mesh->getElementGroup(group).getNodeGroup();
      Array<UInt> new_nodes_list;
      for (const auto & node : nodes_list) {
        if (node_group.find(node) != UInt(-1))
          new_nodes_list.push_back(node);
      }

      this->updateNodalDOFs(dof_id, new_nodes_list);
    }
  }
}

/* -------------------------------------------------------------------------- */
void DOFManager::onNodesRemoved(const Array<UInt> &, const Array<UInt> &,
                                const RemovedNodesEvent &) {}

/* -------------------------------------------------------------------------- */
void DOFManager::onElementsAdded(const Array<Element> &,
                                 const NewElementsEvent &) {}

/* -------------------------------------------------------------------------- */
void DOFManager::onElementsRemoved(const Array<Element> &,
                                   const ElementTypeMapArray<UInt> &,
                                   const RemovedElementsEvent &) {}

/* -------------------------------------------------------------------------- */
void DOFManager::onElementsChanged(const Array<Element> &,
                                   const Array<Element> &,
                                   const ElementTypeMapArray<UInt> &,
                                   const ChangedElementsEvent &) {}

/* -------------------------------------------------------------------------- */

} // namespace akantu
