/**
 * @file   dof_manager.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Aug 12 09:52:30 2015
 *
 * @brief  Implementation of the common parts of the DOFManagers
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "dof_manager.hh"
#include "mesh.hh"
#include "sparse_matrix.hh"
#include "mesh_utils.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
DOFManager::DOFManager(const ID & id, const MemoryID & memory_id)
    : Memory(id, memory_id),
      global_equation_number(0, 1, "global_equation_number"), mesh(NULL),
      local_system_size(0), pure_local_system_size(0), system_size(0) {}

/* -------------------------------------------------------------------------- */
DOFManager::~DOFManager() {
  NodesToElements::scalar_iterator nte_it = this->nodes_to_elements.begin();
  NodesToElements::scalar_iterator nte_end = this->nodes_to_elements.end();
  for (; nte_it != nte_end; ++nte_it)
    delete *nte_it;

  DOFStorage::iterator ds_it = dofs.begin();
  DOFStorage::iterator ds_end = dofs.end();
  for (; ds_it != ds_end; ++ds_it)
    delete ds_it->second;

  SparseMatricesMap::iterator sm_it = matrices.begin();
  SparseMatricesMap::iterator sm_end = matrices.end();
  for (; sm_it != sm_end; ++sm_it)
    delete sm_it->second;

  NonLinearSolversMap::iterator nls_it = non_linear_solvers.begin();
  NonLinearSolversMap::iterator nls_end = non_linear_solvers.end();
  for (; nls_it != nls_end; ++nls_it)
    delete nls_it->second;

  TimeStepSolversMap::iterator tss_it = time_step_solvers.begin();
  TimeStepSolversMap::iterator tss_end = time_step_solvers.end();
  for (; tss_it != tss_end; ++tss_it)
    delete tss_it->second;
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

  UInt * filter_it = NULL;
  if (filter_elements != empty_filter) {
    nb_element = filter_elements.getSize();
    filter_it = filter_elements.storage();
  } else {
    nb_element = this->mesh->getNbElement(type, ghost_type);
  }

  AKANTU_DEBUG_ASSERT(elementary_vect.getSize() == nb_element,
                      "The vector elementary_vect("
                          << elementary_vect.getID()
                          << ") has not the good size.");

  const Array<UInt> connectivity =
      this->mesh->getConnectivity(type, ghost_type);
  Array<UInt>::const_vector_iterator conn_begin =
      connectivity.begin(nb_nodes_per_element);
  Array<UInt>::const_vector_iterator conn_it = conn_begin;

  UInt size_mat = nb_nodes_per_element * nb_degree_of_freedom;
  Array<Real>::const_vector_iterator elem_it = elementary_vect.begin(size_mat);

  for (UInt el = 0; el < nb_element; ++el, ++elem_it) {
    if (filter_it != NULL)
      conn_it = conn_begin + *filter_it;

    for (UInt n = 0, ld = 0; n < nb_nodes_per_element; ++n) {
      UInt offset_node = (*conn_it)(n)*nb_degree_of_freedom;

      for (UInt d = 0; d < nb_degree_of_freedom; ++d, ++ld) {
        array_assembeled[offset_node + d] += scale_factor * (*elem_it)(ld);
      }
    }

    if (filter_it != NULL)
      ++filter_it;
    else
      ++conn_it;
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

  this->assembleToResidual(dof_id, array_localy_assembeled, scale_factor);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DOFManager::registerMesh(Mesh & mesh) {
  this->mesh = &mesh;
  this->mesh->registerEventHandler(*this, 20);

  UInt nb_nodes = this->mesh->getNbNodes();
  this->nodes_to_elements.resize(nb_nodes);
  for (UInt n = 0; n < nb_nodes; ++n) {
    this->nodes_to_elements[n] = new std::set<Element>();
  }
}

/* -------------------------------------------------------------------------- */
DOFManager::DOFData::DOFData() : solution(0, 1, "solution") {}

/* -------------------------------------------------------------------------- */
void DOFManager::registerDOFs(const ID & dof_id, Array<Real> & dofs_array,
                              const DOFSupportType & support_type) {
  DOFStorage::iterator it = this->dofs.find(dof_id);

  if (it != this->dofs.end()) {
    AKANTU_EXCEPTION("This dof array has already been registered");
  }

  DOFData * dofs_storage = new DOFData();
  dofs_storage->dof = &dofs_array;
  dofs_storage->blocked_dofs = NULL;
  dofs_storage->support_type = support_type;

  UInt nb_local_dofs = 0;
  UInt nb_pure_local = 0;

  switch (support_type) {
  case _dst_nodal: {
    nb_local_dofs = this->mesh->getNbNodes();
    AKANTU_DEBUG_ASSERT(
        dofs_array.getSize() == nb_local_dofs,
        "The array of dof is too shot to be associated to nodes.");

    UInt nb_nodes = this->mesh->getNbNodes();
    for (UInt n = 0; n < nb_nodes; ++n) {
      nb_pure_local += this->mesh->isLocalOrMasterNode(n) ? 1 : 0;
    }

    nb_pure_local *= dofs_array.getNbComponent();
    nb_local_dofs *= dofs_array.getNbComponent();
    break;
  }
  case _dst_generic: {
    nb_local_dofs = nb_pure_local =
        dofs_array.getSize() * dofs_array.getNbComponent();
    break;
  }
  default: { AKANTU_EXCEPTION("This type of dofs is not handled yet."); }
  }

  this->pure_local_system_size += nb_pure_local;
  this->local_system_size += nb_local_dofs;

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  comm.allReduce(nb_pure_local, _so_sum);

  this->system_size += nb_pure_local;

  this->dofs[dof_id] = dofs_storage;
}

/* -------------------------------------------------------------------------- */
void DOFManager::registerDOFsDerivative(const ID & dof_id, UInt order,
                                        Array<Real> & dofs_derivative) {
  DOFStorage::iterator it = this->dofs.find(dof_id);

  if (it == this->dofs.end()) {
    AKANTU_EXCEPTION("The dof array corresponding to this derivatives has not "
                     << "been registered yet");
  }

  DOFData & dof = *it->second;
  std::vector<Array<Real> *> & derivatives = dof.dof_derivatives;

  if (derivatives.size() < order) {
    derivatives.resize(order, NULL);
  } else {
    if (derivatives[order - 1] != NULL) {
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
  DOFStorage::iterator it = this->dofs.find(dof_id);

  if (it == this->dofs.end()) {
    AKANTU_EXCEPTION("The dof array corresponding to this derivatives has not "
                     << "been registered yet");
  }

  DOFData & dof = *it->second;

  if (dof.blocked_dofs != NULL) {
    AKANTU_EXCEPTION("The blocked dofs array for "
                     << dof_id << " has already been registered");
  }

  dof.blocked_dofs = &blocked_dofs;
}

/* -------------------------------------------------------------------------- */
void DOFManager::clearJacobian() { this->getMatrix("J").clear(); }

/* -------------------------------------------------------------------------- */
void DOFManager::splitSolutionPerDOFs() {
  DOFStorage::iterator it = this->dofs.begin();
  DOFStorage::iterator end = this->dofs.end();

  for (; it != end; ++it) {
    DOFData & dof_data = *it->second;
    dof_data.solution.resize(dof_data.dof->getSize() *
                             dof_data.dof->getNbComponent());
    this->getSolutionPerDOFs(it->first, dof_data.solution);
  }
}

/* -------------------------------------------------------------------------- */
void DOFManager::registerSparseMatrix(const ID & matrix_id,
                                      SparseMatrix & matrix) {
  SparseMatricesMap::const_iterator it = this->matrices.find(matrix_id);
  if (it != this->matrices.end()) {
    AKANTU_EXCEPTION("The matrix " << matrix_id << " already exists in "
                                   << this->id);
  }

  this->matrices[matrix_id] = &matrix;
}

/* -------------------------------------------------------------------------- */
/// Get an instance of a new SparseMatrix
Array<Real> & DOFManager::getNewLumpedMatrix(const ID & id) {
  ID matrix_id = this->id + ":lumpmtx:" + id;
  LumpedMatricesMap::const_iterator it = this->lumped_matrices.find(matrix_id);
  if (it != this->lumped_matrices.end()) {
    AKANTU_EXCEPTION("The lumped matrix " << matrix_id << " already exists in "
                                          << this->id);
  }

  Array<Real> & mtx = this->alloc<Real>(matrix_id, this->local_system_size, 1);
  this->lumped_matrices[matrix_id] = &mtx;
  return mtx;
}

/* -------------------------------------------------------------------------- */
void DOFManager::registerNonLinearSolver(const ID & non_linear_solver_id,
                                         NonLinearSolver & non_linear_solver) {
  NonLinearSolversMap::const_iterator it =
      this->non_linear_solvers.find(non_linear_solver_id);
  if (it != this->non_linear_solvers.end()) {
    AKANTU_EXCEPTION("The non linear solver " << non_linear_solver_id
                                              << " already exists in "
                                              << this->id);
  }

  this->non_linear_solvers[non_linear_solver_id] = &non_linear_solver;
}

/* -------------------------------------------------------------------------- */
void DOFManager::registerTimeStepSolver(const ID & time_step_solver_id,
                                        TimeStepSolver & time_step_solver) {
  TimeStepSolversMap::const_iterator it =
      this->time_step_solvers.find(time_step_solver_id);
  if (it != this->time_step_solvers.end()) {
    AKANTU_EXCEPTION("The non linear solver " << time_step_solver_id
                                              << " already exists in "
                                              << this->id);
  }

  this->time_step_solvers[time_step_solver_id] = &time_step_solver;
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
Array<Real> & DOFManager::getLumpedMatrix(const ID & id) {
  ID matrix_id = this->id + ":lumpmtx:" + id;
  LumpedMatricesMap::const_iterator it = this->lumped_matrices.find(matrix_id);
  if (it == this->lumped_matrices.end()) {
    AKANTU_SILENT_EXCEPTION("The lumped matrix " << matrix_id << " does not exists in "
                                          << this->id);
  }

  return *(it->second);
}

/* -------------------------------------------------------------------------- */
const Array<Real> & DOFManager::getLumpedMatrix(const ID & id) const {
  ID matrix_id = this->id + ":lumpmtx:" + id;
  LumpedMatricesMap::const_iterator it = this->lumped_matrices.find(matrix_id);
  if (it == this->lumped_matrices.end()) {
    AKANTU_SILENT_EXCEPTION("The lumped matrix " << matrix_id << " does not exists in "
                                          << this->id);
  }

  return *(it->second);
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
void DOFManager::fillNodesToElements() {
  UInt spatial_dimension = this->mesh->getSpatialDimension();
  Element e;

  UInt nb_nodes = this->mesh->getNbNodes();
  for (UInt n = 0; n < nb_nodes; ++n) {
    this->nodes_to_elements[n]->clear();
  }

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {
    Mesh::type_iterator first =
        this->mesh->firstType(spatial_dimension, *gt, _ek_not_defined);
    Mesh::type_iterator last =
        this->mesh->lastType(spatial_dimension, *gt, _ek_not_defined);
    e.ghost_type = *gt;
    for (; first != last; ++first) {
      ElementType type = *first;
      e.type = type;
      e.kind = Mesh::getKind(type);
      UInt nb_element = this->mesh->getNbElement(type, *gt);
      Array<UInt>::const_iterator<Vector<UInt> > conn_it =
          this->mesh->getConnectivity(type, *gt)
              .begin(Mesh::getNbNodesPerElement(type));

      for (UInt el = 0; el < nb_element; ++el, ++conn_it) {
        e.element = el;
        const Vector<UInt> & conn = *conn_it;
        for (UInt n = 0; n < conn.size(); ++n)
          nodes_to_elements[conn(n)]->insert(e);
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
/* Mesh Events                                                                */
/* -------------------------------------------------------------------------- */
void DOFManager::onNodesAdded(const Array<UInt> & nodes_list,
                              __attribute__((unused))
                              const NewNodesEvent & event) {
  Array<UInt>::const_scalar_iterator it = nodes_list.begin();
  Array<UInt>::const_scalar_iterator end = nodes_list.end();

  UInt nb_nodes = this->mesh->getNbNodes();
  this->nodes_to_elements.resize(nb_nodes);

  for (; it != end; ++it) {
    this->nodes_to_elements[*it] = new std::set<Element>();
  }
}

/* -------------------------------------------------------------------------- */
void DOFManager::onNodesRemoved(const Array<UInt> & nodes_list,
                                const Array<UInt> & new_numbering,
                                __attribute__((unused))
                                const RemovedNodesEvent & event) {
  Array<UInt>::const_scalar_iterator it = nodes_list.begin();
  Array<UInt>::const_scalar_iterator end = nodes_list.end();
  for (; it != end; ++it) {
    delete this->nodes_to_elements[*it];
  }

  this->mesh->removeNodesFromArray(this->nodes_to_elements, new_numbering);
}

/* -------------------------------------------------------------------------- */
void DOFManager::onElementsAdded(const Array<Element> & elements_list,
                                 __attribute__((unused))
                                 const NewElementsEvent & event) {
  Array<Element>::const_scalar_iterator it = elements_list.begin();
  Array<Element>::const_scalar_iterator end = elements_list.end();

  for (; it != end; ++it) {
    const Element & elem = *it;
    if (this->mesh->getSpatialDimension(elem.type) !=
        this->mesh->getSpatialDimension())
      continue;

    const Array<UInt> & conn =
        this->mesh->getConnectivity(elem.type, elem.ghost_type);

    UInt nb_nodes_per_elem = this->mesh->getNbNodesPerElement(elem.type);

    for (UInt n = 0; n < nb_nodes_per_elem; ++n) {
      nodes_to_elements[conn(elem.element, n)]->insert(elem);
    }
  }
}

/* -------------------------------------------------------------------------- */
void DOFManager::onElementsRemoved(
    __attribute__((unused)) const Array<Element> & elements_list,
    __attribute__((unused)) const ElementTypeMapArray<UInt> & new_numbering,
    __attribute__((unused)) const RemovedElementsEvent & event) {
  this->fillNodesToElements();
}

/* -------------------------------------------------------------------------- */
void DOFManager::onElementsChanged(
    __attribute__((unused)) const Array<Element> & old_elements_list,
    __attribute__((unused)) const Array<Element> & new_elements_list,
    __attribute__((unused)) const ElementTypeMapArray<UInt> & new_numbering,
    __attribute__((unused)) const ChangedElementsEvent & event) {
  this->fillNodesToElements();
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__
