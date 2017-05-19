/**
 * @file   dof_manager_default.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Tue Aug 11 16:21:01 2015
 *
 * @brief  Implementation of the default DOFManager
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
#include "dof_manager_default.hh"
#include "dof_synchronizer.hh"
#include "element_group.hh"
#include "node_synchronizer.hh"
#include "non_linear_solver_default.hh"
#include "sparse_matrix_aij.hh"
#include "static_communicator.hh"
#include "terms_to_assemble.hh"
#include "time_step_solver_default.hh"
/* -------------------------------------------------------------------------- */
#include <numeric>
#include <unordered_map>
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
inline void DOFManagerDefault::addSymmetricElementalMatrixToSymmetric(
    SparseMatrixAIJ & matrix, const Matrix<Real> & elementary_mat,
    const Vector<UInt> & equation_numbers, UInt max_size) {
  for (UInt i = 0; i < elementary_mat.rows(); ++i) {
    UInt c_irn = equation_numbers(i);
    if (c_irn < max_size) {
      for (UInt j = i; j < elementary_mat.cols(); ++j) {
        UInt c_jcn = equation_numbers(j);
        if (c_jcn < max_size) {
          matrix(c_irn, c_jcn) += elementary_mat(i, j);
        }
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
inline void DOFManagerDefault::addUnsymmetricElementalMatrixToSymmetric(
    SparseMatrixAIJ & matrix, const Matrix<Real> & elementary_mat,
    const Vector<UInt> & equation_numbers, UInt max_size) {
  for (UInt i = 0; i < elementary_mat.rows(); ++i) {
    UInt c_irn = equation_numbers(i);
    if (c_irn < max_size) {
      for (UInt j = 0; j < elementary_mat.cols(); ++j) {
        UInt c_jcn = equation_numbers(j);
        if (c_jcn < max_size) {
          if (c_jcn >= c_irn) {
            matrix(c_irn, c_jcn) += elementary_mat(i, j);
          }
        }
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
inline void DOFManagerDefault::addElementalMatrixToUnsymmetric(
    SparseMatrixAIJ & matrix, const Matrix<Real> & elementary_mat,
    const Vector<UInt> & equation_numbers, UInt max_size) {
  for (UInt i = 0; i < elementary_mat.rows(); ++i) {
    UInt c_irn = equation_numbers(i);
    if (c_irn < max_size) {
      for (UInt j = 0; j < elementary_mat.cols(); ++j) {
        UInt c_jcn = equation_numbers(j);
        if (c_jcn < max_size) {
          matrix(c_irn, c_jcn) += elementary_mat(i, j);
        }
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
DOFManagerDefault::DOFManagerDefault(const ID & id, const MemoryID & memory_id)
    : DOFManager(id, memory_id), residual(0, 1, std::string(id + ":residual")),
      global_residual(nullptr),
      global_solution(0, 1, std::string(id + ":global_solution")),
      global_blocked_dofs(0, 1, std::string(id + ":global_blocked_dofs")),
      previous_global_blocked_dofs(
          0, 1, std::string(id + ":previous_global_blocked_dofs")),
      dofs_type(0, 1, std::string(id + ":dofs_type")),
      data_cache(0, 1, std::string(id + ":data_cache_array")),
      jacobian_release(0),
      global_equation_number(0, 1, "global_equation_number"),
      first_global_dof_id(0), synchronizer(nullptr) {}

/* -------------------------------------------------------------------------- */
DOFManagerDefault::DOFManagerDefault(Mesh & mesh, const ID & id,
                                     const MemoryID & memory_id)
    : DOFManager(mesh, id, memory_id),
      residual(0, 1, std::string(id + ":residual")), global_residual(nullptr),
      global_solution(0, 1, std::string(id + ":global_solution")),
      global_blocked_dofs(0, 1, std::string(id + ":global_blocked_dofs")),
      previous_global_blocked_dofs(
          0, 1, std::string(id + ":previous_global_blocked_dofs")),
      dofs_type(0, 1, std::string(id + ":dofs_type")),
      data_cache(0, 1, std::string(id + ":data_cache_array")),
      jacobian_release(0),
      global_equation_number(0, 1, "global_equation_number"),
      first_global_dof_id(0), synchronizer(nullptr) {
  if (this->mesh->isDistributed())
    this->synchronizer =
        new DOFSynchronizer(*this, this->id + ":dof_synchronizer",
                            this->memory_id, this->mesh->getCommunicator());
}

/* -------------------------------------------------------------------------- */
DOFManagerDefault::~DOFManagerDefault() {
  delete this->synchronizer;
  delete this->global_residual;
}

/* -------------------------------------------------------------------------- */
template <typename T>
void DOFManagerDefault::assembleToGlobalArray(
    const ID & dof_id, const Array<T> & array_to_assemble,
    Array<T> & global_array, T scale_factor) {
  AKANTU_DEBUG_IN();
  const Array<UInt> & equation_number = this->getLocalEquationNumbers(dof_id);

  UInt nb_degree_of_freedoms =
      array_to_assemble.getSize() * array_to_assemble.getNbComponent();

  AKANTU_DEBUG_ASSERT(equation_number.getSize() == nb_degree_of_freedoms,
                      "The array to assemble does not have a correct size."
                          << " (" << array_to_assemble.getID() << ")");

  typename Array<T>::const_scalar_iterator arr_it =
      array_to_assemble.begin_reinterpret(nb_degree_of_freedoms);
  Array<UInt>::const_scalar_iterator equ_it = equation_number.begin();

  for (UInt d = 0; d < nb_degree_of_freedoms; ++d, ++arr_it, ++equ_it) {
    global_array(*equ_it) += scale_factor * (*arr_it);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
DOFManagerDefault::DOFDataDefault::DOFDataDefault(const ID & dof_id)
    : DOFData(dof_id), associated_nodes(0, 1, dof_id + "associated_nodes") {}

/* -------------------------------------------------------------------------- */
DOFManager::DOFData & DOFManagerDefault::getNewDOFData(const ID & dof_id) {
  DOFDataDefault * dofs_storage = new DOFDataDefault(dof_id);
  this->dofs[dof_id] = dofs_storage;
  return *dofs_storage;
}

/* -------------------------------------------------------------------------- */
class GlobalDOFInfoDataAccessor : public DataAccessor<UInt> {
public:
  typedef
      typename std::unordered_map<UInt, std::vector<UInt>>::size_type size_type;

  GlobalDOFInfoDataAccessor() = default;

  void addDOFToNode(UInt node, UInt dof) { dofs_per_node[node].push_back(dof); }
  UInt getNthDOFForNode(UInt nth_dof, UInt node) const {
    return dofs_per_node.find(node)->second[nth_dof];
  }

  virtual UInt getNbData(const Array<UInt> & nodes,
                         const SynchronizationTag & tag) const {
    if (tag == _gst_size) {
      return nodes.getSize() * sizeof(size_type);
    }

    if (tag == _gst_update) {
      UInt total_size = 0;
      for (auto node : nodes) {
        auto it = dofs_per_node.find(node);
        if (it != dofs_per_node.end())
          total_size += CommunicationBuffer::sizeInBuffer(it->second);
      }
      return total_size;
    }

    return 0;
  }

  virtual void packData(CommunicationBuffer & buffer, const Array<UInt> & nodes,
                        const SynchronizationTag & tag) const {
    for (auto node : nodes) {
      auto it = dofs_per_node.find(node);
      if (tag == _gst_size) {
        if (it != dofs_per_node.end()) {
          buffer << it->second.size();
        } else {
          buffer << 0;
        }
      } else if (tag == _gst_update) {
        if (it != dofs_per_node.end())
          buffer << it->second;
      }
    }
  }

  virtual void unpackData(CommunicationBuffer & buffer,
                          const Array<UInt> & nodes,
                          const SynchronizationTag & tag) {
    for (auto node : nodes) {
      auto it = dofs_per_node.find(node);
      if (tag == _gst_size) {
        size_type size;
        buffer >> size;
        if (size != 0)
          dofs_per_node[node].resize(size);
      } else if (tag == _gst_update) {
        if (it != dofs_per_node.end())
          buffer >> it->second;
      }
    }
  }

protected:
  std::unordered_map<UInt, std::vector<UInt>> dofs_per_node;
};

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::registerDOFsInternal(const ID & dof_id, UInt nb_dofs,
                                             UInt nb_pure_local_dofs) {
  // Count the number of pure local dofs per proc
  const StaticCommunicator * comm = nullptr;
  UInt prank = 0;
  UInt psize = 1;

  if (mesh) {
    comm = &mesh->getCommunicator();
    prank = comm->whoAmI();
    psize = comm->getNbProc();
  } else {
    comm = &(StaticCommunicator::getStaticCommunicator());
  }

  // access the relevant data to update
  DOFDataDefault & dof_data = this->getDOFDataTyped<DOFDataDefault>(dof_id);
  const DOFSupportType & support_type = dof_data.support_type;
  const ID & group = dof_data.group_support;

  GlobalDOFInfoDataAccessor data_accessor;

  // resize all relevant arrays
  this->residual.resize(this->local_system_size);
  this->dofs_type.resize(local_system_size);
  this->global_solution.resize(this->local_system_size);
  this->global_blocked_dofs.resize(this->local_system_size);
  this->previous_global_blocked_dofs.resize(this->local_system_size);
  this->global_equation_number.resize(this->local_system_size);
  dof_data.local_equation_number.resize(nb_dofs);

  // determine the first local/global dof id to use
  Array<UInt> nb_dofs_per_proc(psize);
  nb_dofs_per_proc(prank) = nb_pure_local_dofs;
  comm->allGather(nb_dofs_per_proc);
  this->first_global_dof_id += std::accumulate(
      nb_dofs_per_proc.begin(), nb_dofs_per_proc.begin() + prank, 0);
  UInt first_dof_id = this->local_system_size - nb_dofs;

  const Array<UInt> * support_nodes = nullptr;
  if (support_type == _dst_nodal) {
    if (group != "mesh") {
      support_nodes =
          &this->mesh->getElementGroup(group).getNodeGroup().getNodes();
    }

    dof_data.associated_nodes.resize(nb_dofs);
  }

  // update per dof info
  for (UInt d = 0; d < nb_dofs; ++d) {
    UInt local_eq_num = first_dof_id + d;
    dof_data.local_equation_number(d) = local_eq_num;

    bool is_local_dof = true;

    // determine the dof type
    switch (support_type) {
    case _dst_nodal: {
      UInt node = d / dof_data.dof->getNbComponent();

      if (support_nodes)
        node = (*support_nodes)(node);

      this->dofs_type(local_eq_num) = this->mesh->getNodeType(node);
      dof_data.associated_nodes(d) = node;

      is_local_dof = this->mesh->isLocalOrMasterNode(node);

      if (is_local_dof) {
        data_accessor.addDOFToNode(node, first_global_dof_id);
      }
      break;
    }
    case _dst_generic: {
      this->dofs_type(local_eq_num) = _nt_normal;
      break;
    }
    default: { AKANTU_EXCEPTION("This type of dofs is not handled yet."); }
    }

    // update global id for local dofs
    if (is_local_dof) {
      this->global_equation_number(local_eq_num) = this->first_global_dof_id;
      this->global_to_local_mapping[this->first_global_dof_id] = local_eq_num;
      ++this->first_global_dof_id;
    } else {
      this->global_equation_number(local_eq_num) = 0;
    }
  }

  if (support_type == _dst_nodal) {
    auto & node_synchronizer = this->mesh->getNodeSynchronizer();
    node_synchronizer.synchronizeOnce(data_accessor, _gst_size);
    node_synchronizer.synchronizeOnce(data_accessor, _gst_update);

    std::vector<UInt> counters(nb_dofs);

    for (UInt d = 0; d < nb_dofs; ++d) {
      UInt local_eq_num = first_dof_id + d;
      if (this->isSlaveDOF(local_eq_num)) {
        UInt node = d / dof_data.dof->getNbComponent();

        UInt dof_count = counters[node]++;
        UInt global_dof_id = data_accessor.getNthDOFForNode(dof_count, node);

        this->global_equation_number(local_eq_num) = global_dof_id;
        this->global_to_local_mapping[global_dof_id] = local_eq_num;
      }
    }
  }
  // update the synchronizer if needed
  if (this->synchronizer)
    this->synchronizer->registerDOFs(dof_id);
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::registerDOFs(const ID & dof_id,
                                     Array<Real> & dofs_array,
                                     const DOFSupportType & support_type) {
  // stores the current numbers of dofs
  UInt nb_dofs_old = this->local_system_size;
  UInt nb_pure_local_dofs_old = this->pure_local_system_size;

  // update or create the dof_data
  DOFManager::registerDOFs(dof_id, dofs_array, support_type);

  UInt nb_dofs = this->local_system_size - nb_dofs_old;
  UInt nb_pure_local_dofs =
      this->pure_local_system_size - nb_pure_local_dofs_old;

  this->registerDOFsInternal(dof_id, nb_dofs, nb_pure_local_dofs);
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::registerDOFs(const ID & dof_id,
                                     Array<Real> & dofs_array,
                                     const ID & group_support) {
  // stores the current numbers of dofs
  UInt nb_dofs_old = this->local_system_size;
  UInt nb_pure_local_dofs_old = this->pure_local_system_size;

  // update or create the dof_data
  DOFManager::registerDOFs(dof_id, dofs_array, group_support);

  UInt nb_dofs = this->local_system_size - nb_dofs_old;
  UInt nb_pure_local_dofs =
      this->pure_local_system_size - nb_pure_local_dofs_old;

  this->registerDOFsInternal(dof_id, nb_dofs, nb_pure_local_dofs);
}

/* -------------------------------------------------------------------------- */
SparseMatrix & DOFManagerDefault::getNewMatrix(const ID & id,
                                               const MatrixType & matrix_type) {
  ID matrix_id = this->id + ":mtx:" + id;
  SparseMatrix * sm = new SparseMatrixAIJ(*this, matrix_type, matrix_id);
  this->registerSparseMatrix(matrix_id, *sm);

  return *sm;
}

/* -------------------------------------------------------------------------- */
SparseMatrix & DOFManagerDefault::getNewMatrix(const ID & id,
                                               const ID & matrix_to_copy_id) {

  ID matrix_id = this->id + ":mtx:" + id;
  SparseMatrixAIJ & sm_to_copy = this->getMatrix(matrix_to_copy_id);
  SparseMatrix * sm = new SparseMatrixAIJ(sm_to_copy, matrix_id);
  this->registerSparseMatrix(matrix_id, *sm);

  return *sm;
}

/* -------------------------------------------------------------------------- */
SparseMatrixAIJ & DOFManagerDefault::getMatrix(const ID & id) {
  SparseMatrix & matrix = DOFManager::getMatrix(id);

  return dynamic_cast<SparseMatrixAIJ &>(matrix);
}

/* -------------------------------------------------------------------------- */
NonLinearSolver &
DOFManagerDefault::getNewNonLinearSolver(const ID & id,
                                         const NonLinearSolverType & type) {
  ID non_linear_solver_id = this->id + ":nls:" + id;
  NonLinearSolver * nls = NULL;
  switch (type) {
#if defined(AKANTU_IMPLICIT)
  case _nls_newton_raphson:
  case _nls_newton_raphson_modified: {
    nls = new NonLinearSolverNewtonRaphson(*this, type, non_linear_solver_id,
                                           this->memory_id);
    break;
  }
  case _nls_linear: {
    nls = new NonLinearSolverLinear(*this, type, non_linear_solver_id,
                                    this->memory_id);
    break;
  }
#endif
  case _nls_lumped: {
    nls = new NonLinearSolverLumped(*this, type, non_linear_solver_id,
                                    this->memory_id);
    break;
  }
  default:
    AKANTU_EXCEPTION("The asked type of non linear solver is not supported by "
                     "this dof manager");
  }

  this->registerNonLinearSolver(non_linear_solver_id, *nls);

  return *nls;
}

/* -------------------------------------------------------------------------- */
TimeStepSolver &
DOFManagerDefault::getNewTimeStepSolver(const ID & id,
                                        const TimeStepSolverType & type,
                                        NonLinearSolver & non_linear_solver) {
  ID time_step_solver_id = this->id + ":tss:" + id;

  TimeStepSolver * tss = new TimeStepSolverDefault(
      *this, type, non_linear_solver, time_step_solver_id, this->memory_id);

  this->registerTimeStepSolver(time_step_solver_id, *tss);

  return *tss;
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::clearResidual() {
  this->residual.resize(this->local_system_size);
  this->residual.clear();
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::clearMatrix(const ID & mtx) {
  this->getMatrix(mtx).clear();
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::clearLumpedMatrix(const ID & mtx) {
  this->getLumpedMatrix(mtx).clear();
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::updateGlobalBlockedDofs() {
  DOFStorage::iterator it = this->dofs.begin();
  DOFStorage::iterator end = this->dofs.end();

  this->previous_global_blocked_dofs.copy(this->global_blocked_dofs);
  this->global_blocked_dofs.resize(this->local_system_size);
  this->global_blocked_dofs.clear();

  for (; it != end; ++it) {
    if (!this->hasBlockedDOFs(it->first)) continue;

    DOFData & dof_data = *it->second;
    this->assembleToGlobalArray(it->first, *dof_data.blocked_dofs,
                                this->global_blocked_dofs, true);
  }
}

/* -------------------------------------------------------------------------- */
template <typename T>
void DOFManagerDefault::getArrayPerDOFs(const ID & dof_id,
                                        const Array<T> & global_array,
                                        Array<T> & local_array) const {
  AKANTU_DEBUG_IN();

  const Array<UInt> & equation_number = this->getLocalEquationNumbers(dof_id);

  UInt nb_degree_of_freedoms = equation_number.getSize();
  local_array.resize(nb_degree_of_freedoms / local_array.getNbComponent());

  auto loc_it = local_array.begin_reinterpret(nb_degree_of_freedoms);
  auto equ_it = equation_number.begin();

  for (UInt d = 0; d < nb_degree_of_freedoms; ++d, ++loc_it, ++equ_it) {
    (*loc_it) = global_array(*equ_it);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::getEquationsNumbers(const ID & dof_id,
                                            Array<UInt> & equation_numbers) {
  AKANTU_DEBUG_IN();
  this->getArrayPerDOFs(dof_id, this->global_equation_number, equation_numbers);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::getSolutionPerDOFs(const ID & dof_id,
                                           Array<Real> & solution_array) {
  AKANTU_DEBUG_IN();
  this->getArrayPerDOFs(dof_id, this->global_solution, solution_array);
  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
void DOFManagerDefault::getLumpedMatrixPerDOFs(const ID & dof_id,
                                               const ID & lumped_mtx,
                                               Array<Real> & lumped) {
  AKANTU_DEBUG_IN();
  this->getArrayPerDOFs(dof_id, this->getLumpedMatrix(lumped_mtx), lumped);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::assembleToResidual(
    const ID & dof_id, const Array<Real> & array_to_assemble,
    Real scale_factor) {
  AKANTU_DEBUG_IN();

  this->assembleToGlobalArray(dof_id, array_to_assemble, this->residual,
                              scale_factor);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::assembleToLumpedMatrix(
    const ID & dof_id, const Array<Real> & array_to_assemble,
    const ID & lumped_mtx, Real scale_factor) {
  AKANTU_DEBUG_IN();

  Array<Real> & lumped = this->getLumpedMatrix(lumped_mtx);
  this->assembleToGlobalArray(dof_id, array_to_assemble, lumped, scale_factor);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::assembleMatMulVectToResidual(const ID & dof_id,
                                                     const ID & A_id,
                                                     const Array<Real> & x,
                                                     Real scale_factor) {
  SparseMatrixAIJ & A = this->getMatrix(A_id);

  // Array<Real> data_cache(this->local_system_size, 1, 0.);
  this->data_cache.resize(this->local_system_size);
  this->data_cache.clear();

  this->assembleToGlobalArray(dof_id, x, data_cache, 1.);

  A.matVecMul(data_cache, this->residual, scale_factor, 1.);
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::assembleLumpedMatMulVectToResidual(
    const ID & dof_id, const ID & A_id, const Array<Real> & x,
    Real scale_factor) {
  const Array<Real> & A = this->getLumpedMatrix(A_id);

  //  Array<Real> data_cache(this->local_system_size, 1, 0.);
  this->data_cache.resize(this->local_system_size);
  this->data_cache.clear();

  this->assembleToGlobalArray(dof_id, x, data_cache, scale_factor);

  Array<Real>::const_scalar_iterator A_it = A.begin();
  Array<Real>::const_scalar_iterator A_end = A.end();
  Array<Real>::const_scalar_iterator x_it = data_cache.begin();
  Array<Real>::scalar_iterator r_it = this->residual.begin();

  for (; A_it != A_end; ++A_it, ++x_it, ++r_it) {
    *r_it += *A_it * *x_it;
  }
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::assembleElementalMatricesToMatrix(
    const ID & matrix_id, const ID & dof_id, const Array<Real> & elementary_mat,
    const ElementType & type, const GhostType & ghost_type,
    const MatrixType & elemental_matrix_type,
    const Array<UInt> & filter_elements) {
  AKANTU_DEBUG_IN();

  this->addToProfile(matrix_id, dof_id, type, ghost_type);

  DOFData & dof_data = this->getDOFData(dof_id);

  const Array<UInt> & equation_number = this->getLocalEquationNumbers(dof_id);
  SparseMatrixAIJ & A = this->getMatrix(matrix_id);

  UInt nb_element;
  if (ghost_type == _not_ghost) {
    nb_element = this->mesh->getNbElement(type);
  } else {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  UInt * filter_it = nullptr;
  if (filter_elements != empty_filter) {
    nb_element = filter_elements.getSize();
    filter_it = filter_elements.storage();
  } else {
    if (dof_data.group_support != "mesh") {
      const Array<UInt> & group_elements =
          this->mesh->getElementGroup(dof_data.group_support)
              .getElements(type, ghost_type);
      nb_element = group_elements.getSize();
      filter_it = group_elements.storage();
    } else {
      nb_element = this->mesh->getNbElement(type, ghost_type);
    }
  }

  AKANTU_DEBUG_ASSERT(elementary_mat.getSize() == nb_element,
                      "The vector elementary_mat("
                          << elementary_mat.getID()
                          << ") has not the good size.");

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

  UInt nb_degree_of_freedom = dof_data.dof->getNbComponent();

  const Array<UInt> & connectivity =
      this->mesh->getConnectivity(type, ghost_type);
  auto conn_begin = connectivity.begin(nb_nodes_per_element);
  auto conn_it = conn_begin;

  UInt size_mat = nb_nodes_per_element * nb_degree_of_freedom;

  Vector<UInt> element_eq_nb(nb_degree_of_freedom * nb_nodes_per_element);
  Array<Real>::const_matrix_iterator el_mat_it =
      elementary_mat.begin(size_mat, size_mat);

  for (UInt e = 0; e < nb_element; ++e, ++el_mat_it) {
    if (filter_it)
      conn_it = conn_begin + *filter_it;

    this->extractElementEquationNumber(equation_number, *conn_it,
                                       nb_degree_of_freedom, element_eq_nb);
    std::transform(element_eq_nb.storage(),
                   element_eq_nb.storage() + element_eq_nb.size(),
                   element_eq_nb.storage(), [&](UInt & local) -> UInt {
                     return this->localToGlobalEquationNumber(local);
                   });

    if (filter_it)
      ++filter_it;
    else
      ++conn_it;

    if (A.getMatrixType() == _symmetric)
      if (elemental_matrix_type == _symmetric)
        this->addSymmetricElementalMatrixToSymmetric(
            A, *el_mat_it, element_eq_nb, A.getSize());
      else
        this->addUnsymmetricElementalMatrixToSymmetric(
            A, *el_mat_it, element_eq_nb, A.getSize());
    else
      this->addElementalMatrixToUnsymmetric(A, *el_mat_it, element_eq_nb,
                                            A.getSize());
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::assemblePreassembledMatrix(
    const ID & dof_id_m, const ID & dof_id_n, const ID & matrix_id,
    const TermsToAssemble & terms) {
  const Array<UInt> & equation_number_m =
      this->getLocalEquationNumbers(dof_id_m);
  const Array<UInt> & equation_number_n =
      this->getLocalEquationNumbers(dof_id_n);
  SparseMatrixAIJ & A = this->getMatrix(matrix_id);

  for (const auto & term : terms) {
    A.addToMatrix(equation_number_m(term.i()), equation_number_n(term.j()),
                  term);
  }
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::addToProfile(const ID & matrix_id, const ID & dof_id,
                                     const ElementType & type,
                                     const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  const DOFData & dof_data = this->getDOFData(dof_id);

  if (dof_data.support_type != _dst_nodal)
    return;

  auto mat_dof = std::make_pair(matrix_id, dof_id);
  auto type_pair = std::make_pair(type, ghost_type);

  auto prof_it = this->matrix_profiled_dofs.find(mat_dof);
  if (prof_it != this->matrix_profiled_dofs.end() &&
      std::find(prof_it->second.begin(), prof_it->second.end(), type_pair) !=
          prof_it->second.end())
    return;

  UInt nb_degree_of_freedom_per_node = dof_data.dof->getNbComponent();

  const Array<UInt> & equation_number = this->getLocalEquationNumbers(dof_id);

  SparseMatrixAIJ & A = this->getMatrix(matrix_id);

  UInt size = A.getSize();

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

  const auto & connectivity = this->mesh->getConnectivity(type, ghost_type);
  auto cbegin = connectivity.begin(nb_nodes_per_element);
  auto cit = cbegin;

  UInt nb_elements = connectivity.getSize();
  UInt * ge_it = nullptr;
  if (dof_data.group_support != "mesh") {
    const Array<UInt> & group_elements =
        this->mesh->getElementGroup(dof_data.group_support)
            .getElements(type, ghost_type);
    ge_it = group_elements.storage();
    nb_elements = group_elements.getSize();
  }

  UInt size_mat = nb_nodes_per_element * nb_degree_of_freedom_per_node;
  Vector<UInt> element_eq_nb(size_mat);

  for (UInt e = 0; e < nb_elements; ++e) {
    if (ge_it)
      cit = cbegin + *ge_it;

    this->extractElementEquationNumber(
        equation_number, *cit, nb_degree_of_freedom_per_node, element_eq_nb);
    std::transform(element_eq_nb.storage(),
                   element_eq_nb.storage() + element_eq_nb.size(),
                   element_eq_nb.storage(), [&](UInt & local) -> UInt {
                     return this->localToGlobalEquationNumber(local);
                   });

    if (ge_it)
      ++ge_it;
    else
      ++cit;

    for (UInt i = 0; i < size_mat; ++i) {
      UInt c_irn = element_eq_nb(i);
      if (c_irn < size) {
        for (UInt j = 0; j < size_mat; ++j) {
          UInt c_jcn = element_eq_nb(j);
          if (c_jcn < size) {
            A.addToProfile(c_irn, c_jcn);
          }
        }
      }
    }
  }

  this->matrix_profiled_dofs[mat_dof].push_back(type_pair);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::applyBoundary(const ID & matrix_id) {
  this->updateGlobalBlockedDofs();
  SparseMatrixAIJ & J = this->getMatrix(matrix_id);

  if (this->jacobian_release == J.getValueRelease()) {
    Array<bool>::const_scalar_iterator it = global_blocked_dofs.begin();
    Array<bool>::const_scalar_iterator end = global_blocked_dofs.end();

    Array<bool>::const_scalar_iterator pit =
        previous_global_blocked_dofs.begin();

    for (; it != end && *it == *pit; ++it, ++pit)
      ;

    if (it != end)
      J.applyBoundary();
  } else {
    J.applyBoundary();
  }

  this->jacobian_release = J.getValueRelease();
}

/* -------------------------------------------------------------------------- */
const Array<Real> & DOFManagerDefault::getGlobalResidual() {
  if (this->synchronizer) {
    if (not this->global_residual) {
      this->global_residual = new Array<Real>(0, 1, "global_residual");
    }

    if (this->synchronizer->getCommunicator().whoAmI() == 0) {
      this->global_residual->resize(this->system_size);
      this->synchronizer->gather(this->residual, *this->global_residual);
    } else {
      this->synchronizer->gather(this->residual);
    }

    return *this->global_residual;
  } else {
    return this->residual;
  }
}

/* -------------------------------------------------------------------------- */
const Array<Real> & DOFManagerDefault::getResidual() const {
  return this->residual;
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::setGlobalSolution(const Array<Real> & solution) {
  if (this->synchronizer) {
    if (this->synchronizer->getCommunicator().whoAmI() == 0) {
      this->synchronizer->scatter(this->global_solution, solution);
    } else {
      this->synchronizer->scatter(this->global_solution);
    }
  } else {
    AKANTU_DEBUG_ASSERT(solution.getSize() == this->global_solution.getSize(),
                        "Sequential call to this function needs the solution "
                        "to be the same size as the global_solution");
    this->global_solution.copy(solution);
  }
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__

//  LocalWords:  dof dofs
