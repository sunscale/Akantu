/**
 * @file   dof_manager_default.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Aug 18 2015
 * @date last modification: Thu Feb 08 2018
 *
 * @brief  Implementation of the default DOFManager
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
#include "dof_manager_default.hh"
#include "communicator.hh"
#include "dof_synchronizer.hh"
#include "element_group.hh"
#include "node_synchronizer.hh"
#include "non_linear_solver_default.hh"
#include "sparse_matrix_aij.hh"
#include "terms_to_assemble.hh"
#include "time_step_solver_default.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
#include <memory>
#include <numeric>
#include <unordered_map>
/* -------------------------------------------------------------------------- */

namespace akantu {

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
      dofs_flag(0, 1, std::string(id + ":dofs_type")),
      data_cache(0, 1, std::string(id + ":data_cache_array")),
      global_equation_number(0, 1, "global_equation_number"),
      synchronizer(nullptr) {}

/* -------------------------------------------------------------------------- */
DOFManagerDefault::DOFManagerDefault(Mesh & mesh, const ID & id,
                                     const MemoryID & memory_id)
    : DOFManager(mesh, id, memory_id),
      residual(0, 1, std::string(id + ":residual")), global_residual(nullptr),
      global_solution(0, 1, std::string(id + ":global_solution")),
      global_blocked_dofs(0, 1, std::string(id + ":global_blocked_dofs")),
      previous_global_blocked_dofs(
          0, 1, std::string(id + ":previous_global_blocked_dofs")),
      dofs_flag(0, 1, std::string(id + ":dofs_type")),
      data_cache(0, 1, std::string(id + ":data_cache_array")),
      global_equation_number(0, 1, "global_equation_number"),
      first_global_dof_id(0), synchronizer(nullptr) {
  if (this->mesh->isDistributed())
    this->synchronizer = std::make_unique<DOFSynchronizer>(
        *this, this->id + ":dof_synchronizer", this->memory_id);
}

/* -------------------------------------------------------------------------- */
DOFManagerDefault::~DOFManagerDefault() = default;

/* -------------------------------------------------------------------------- */
template <typename T>
void DOFManagerDefault::assembleToGlobalArray(
    const ID & dof_id, const Array<T> & array_to_assemble,
    Array<T> & global_array, T scale_factor) {
  AKANTU_DEBUG_IN();
  const Array<UInt> & equation_number = this->getLocalEquationNumbers(dof_id);

  UInt nb_degree_of_freedoms =
      array_to_assemble.size() * array_to_assemble.getNbComponent();

  AKANTU_DEBUG_ASSERT(equation_number.size() == nb_degree_of_freedoms,
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
  this->dofs[dof_id] = std::make_unique<DOFDataDefault>(dof_id);
  return *this->dofs[dof_id];
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::registerDOFsInternal(const ID & dof_id, UInt nb_dofs,
                                             UInt nb_pure_local_dofs) {
  // auto prank = this->communicator.whoAmI();
  // auto psize = this->communicator.getNbProc();

  // access the relevant data to update
  auto & dof_data = this->getDOFDataTyped<DOFDataDefault>(dof_id);
  const auto & support_type = dof_data.support_type;

  const auto & group = dof_data.group_support;

  if (support_type == _dst_nodal and group != "__mesh__") {
    auto & support_nodes =
        this->mesh->getElementGroup(group).getNodeGroup().getNodes();
    this->updateDOFsData(
        dof_data, nb_dofs, nb_pure_local_dofs,
        [&support_nodes](UInt node) -> UInt { return support_nodes[node]; });
  } else {

    this->updateDOFsData(dof_data, nb_dofs, nb_pure_local_dofs,
                         [](UInt node) -> UInt { return node; });
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
  std::unique_ptr<SparseMatrix> sm =
      std::make_unique<SparseMatrixAIJ>(*this, matrix_type, matrix_id);
  return this->registerSparseMatrix(matrix_id, sm);
}

/* -------------------------------------------------------------------------- */
SparseMatrix & DOFManagerDefault::getNewMatrix(const ID & id,
                                               const ID & matrix_to_copy_id) {

  ID matrix_id = this->id + ":mtx:" + id;
  SparseMatrixAIJ & sm_to_copy = this->getMatrix(matrix_to_copy_id);
  std::unique_ptr<SparseMatrix> sm =
      std::make_unique<SparseMatrixAIJ>(sm_to_copy, matrix_id);
  return this->registerSparseMatrix(matrix_id, sm);
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
  std::unique_ptr<NonLinearSolver> nls;
  switch (type) {
#if defined(AKANTU_IMPLICIT)
  case _nls_newton_raphson:
  case _nls_newton_raphson_modified: {
    nls = std::make_unique<NonLinearSolverNewtonRaphson>(
        *this, type, non_linear_solver_id, this->memory_id);
    break;
  }
  case _nls_linear: {
    nls = std::make_unique<NonLinearSolverLinear>(
        *this, type, non_linear_solver_id, this->memory_id);
    break;
  }
#endif
  case _nls_lumped: {
    nls = std::make_unique<NonLinearSolverLumped>(
        *this, type, non_linear_solver_id, this->memory_id);
    break;
  }
  default:
    AKANTU_EXCEPTION("The asked type of non linear solver is not supported by "
                     "this dof manager");
  }

  return this->registerNonLinearSolver(non_linear_solver_id, nls);
}

/* -------------------------------------------------------------------------- */
TimeStepSolver &
DOFManagerDefault::getNewTimeStepSolver(const ID & id,
                                        const TimeStepSolverType & type,
                                        NonLinearSolver & non_linear_solver) {
  ID time_step_solver_id = this->id + ":tss:" + id;

  std::unique_ptr<TimeStepSolver> tss = std::make_unique<TimeStepSolverDefault>(
      *this, type, non_linear_solver, time_step_solver_id, this->memory_id);

  return this->registerTimeStepSolver(time_step_solver_id, tss);
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
  auto it = this->dofs.begin();
  auto end = this->dofs.end();

  this->previous_global_blocked_dofs.copy(this->global_blocked_dofs);
  this->global_blocked_dofs.resize(this->local_system_size);
  this->global_blocked_dofs.clear();

  for (; it != end; ++it) {
    if (!this->hasBlockedDOFs(it->first))
      continue;

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

  UInt nb_degree_of_freedoms = equation_number.size();
  local_array.resize(nb_degree_of_freedoms / local_array.getNbComponent());

  auto loc_it = local_array.begin_reinterpret(nb_degree_of_freedoms);
  auto equ_it = equation_number.begin();

  for (UInt d = 0; d < nb_degree_of_freedoms; ++d, ++loc_it, ++equ_it) {
    (*loc_it) = global_array(*equ_it);
  }

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

  Array<Real> tmp_residual(this->residual.size(), 1, 0.);

  A.matVecMul(data_cache, tmp_residual, scale_factor, 1.);
  this->residual += tmp_residual;
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

  auto A_it = A.begin();
  auto A_end = A.end();
  auto x_it = data_cache.begin();
  auto r_it = this->residual.begin();

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

  DOFData & dof_data = this->getDOFData(dof_id);

  AKANTU_DEBUG_ASSERT(dof_data.support_type == _dst_nodal,
                      "This function applies only on Nodal dofs");

  this->addToProfile(matrix_id, dof_id, type, ghost_type);

  const Array<UInt> & equation_number = this->getLocalEquationNumbers(dof_id);
  SparseMatrixAIJ & A = this->getMatrix(matrix_id);

  UInt nb_element;
  UInt * filter_it = nullptr;
  if (filter_elements != empty_filter) {
    nb_element = filter_elements.size();
    filter_it = filter_elements.storage();
  } else {
    if (dof_data.group_support != "__mesh__") {
      const Array<UInt> & group_elements =
          this->mesh->getElementGroup(dof_data.group_support)
              .getElements(type, ghost_type);
      nb_element = group_elements.size();
      filter_it = group_elements.storage();
    } else {
      nb_element = this->mesh->getNbElement(type, ghost_type);
    }
  }

  AKANTU_DEBUG_ASSERT(elementary_mat.size() == nb_element,
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
    std::transform(element_eq_nb.begin(), element_eq_nb.end(),
                   element_eq_nb.begin(), [&](UInt & local) -> UInt {
                     return this->localToGlobalEquationNumber(local);
                   });

    if (filter_it)
      ++filter_it;
    else
      ++conn_it;

    if (A.getMatrixType() == _symmetric)
      if (elemental_matrix_type == _symmetric)
        this->addSymmetricElementalMatrixToSymmetric(A, *el_mat_it,
                                                     element_eq_nb, A.size());
      else
        this->addUnsymmetricElementalMatrixToSymmetric(A, *el_mat_it,
                                                       element_eq_nb, A.size());
    else
      this->addElementalMatrixToUnsymmetric(A, *el_mat_it, element_eq_nb,
                                            A.size());
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
    UInt gi = this->localToGlobalEquationNumber(equation_number_m(term.i()));
    UInt gj = this->localToGlobalEquationNumber(equation_number_n(term.j()));
    A.add(gi, gj, term);
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

  UInt size = A.size();

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

  const auto & connectivity = this->mesh->getConnectivity(type, ghost_type);
  auto cbegin = connectivity.begin(nb_nodes_per_element);
  auto cit = cbegin;

  UInt nb_elements = connectivity.size();
  UInt * ge_it = nullptr;
  if (dof_data.group_support != "__mesh__") {
    const Array<UInt> & group_elements =
        this->mesh->getElementGroup(dof_data.group_support)
            .getElements(type, ghost_type);
    ge_it = group_elements.storage();
    nb_elements = group_elements.size();
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
            A.add(c_irn, c_jcn);
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
  SparseMatrixAIJ & J = this->getMatrix(matrix_id);

  if (this->jacobian_release == J.getValueRelease()) {
    auto are_equal =
        std::equal(global_blocked_dofs.begin(), global_blocked_dofs.end(),
                   previous_global_blocked_dofs.begin());

    if (not are_equal)
      J.applyBoundary();

    previous_global_blocked_dofs.copy(global_blocked_dofs);
  } else {
    J.applyBoundary();
  }

  this->jacobian_release = J.getValueRelease();
}

/* -------------------------------------------------------------------------- */
const Array<Real> & DOFManagerDefault::getGlobalResidual() {
  if (this->synchronizer) {
    if (not this->global_residual) {
      this->global_residual =
          std::make_unique<Array<Real>>(0, 1, "global_residual");
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
    AKANTU_DEBUG_ASSERT(solution.size() == this->global_solution.size(),
                        "Sequential call to this function needs the solution "
                        "to be the same size as the global_solution");
    this->global_solution.copy(solution);
  }
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::onNodesAdded(const Array<UInt> & nodes_list,
                                     const NewNodesEvent & event) {
  DOFManager::onNodesAdded(nodes_list, event);

  if (this->synchronizer)
    this->synchronizer->onNodesAdded(nodes_list);
}

/* -------------------------------------------------------------------------- */
std::pair<UInt, UInt>
DOFManagerDefault::updateNodalDOFs(const ID & dof_id,
                                   const Array<UInt> & nodes_list) {
  UInt nb_new_local_dofs, nb_new_pure_local;
  std::tie(nb_new_local_dofs, nb_new_pure_local) =
      DOFManager::updateNodalDOFs(dof_id, nodes_list);

  auto & dof_data = this->getDOFDataTyped<DOFDataDefault>(dof_id);
  updateDOFsData(dof_data, nb_new_local_dofs, nb_new_pure_local,
                 [&nodes_list](UInt pos) -> UInt { return nodes_list[pos]; });

  return std::make_pair(nb_new_local_dofs, nb_new_pure_local);
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
class GlobalDOFInfoDataAccessor : public DataAccessor<UInt> {
public:
  using size_type =
      typename std::unordered_map<UInt, std::vector<UInt>>::size_type;

  GlobalDOFInfoDataAccessor(DOFManagerDefault::DOFDataDefault & dof_data,
                            DOFManagerDefault & dof_manager)
      : dof_data(dof_data), dof_manager(dof_manager) {
    for (auto && pair :
         zip(dof_data.local_equation_number, dof_data.associated_nodes)) {
      UInt node, dof;
      std::tie(dof, node) = pair;
      dofs_per_node[node].push_back(dof);
    }
  }

  UInt getNbData(const Array<UInt> & nodes,
                 const SynchronizationTag & tag) const override {
    if (tag == _gst_ask_nodes) {
      return nodes.size() * dof_data.dof->getNbComponent() * sizeof(UInt);
    }

    return 0;
  }

  void packData(CommunicationBuffer & buffer, const Array<UInt> & nodes,
                const SynchronizationTag & tag) const override {
    if (tag == _gst_ask_nodes) {
      for (auto & node : nodes) {
        auto & dofs = dofs_per_node.at(node);
        for (auto & dof : dofs) {
          buffer << dof_manager.global_equation_number(dof);
        }
      }
    }
  }

  void unpackData(CommunicationBuffer & buffer, const Array<UInt> & nodes,
                  const SynchronizationTag & tag) override {
    if (tag == _gst_ask_nodes) {
      for (auto & node : nodes) {
        auto & dofs = dofs_per_node[node];
        for (auto dof : dofs) {
          UInt global_dof;
          buffer >> global_dof;
          dof_manager.global_equation_number(dof) = global_dof;
          dof_manager.global_to_local_mapping[global_dof] = dof;
        }
      }
    }
  }

protected:
  std::unordered_map<UInt, std::vector<UInt>> dofs_per_node;
  DOFManagerDefault::DOFDataDefault & dof_data;
  DOFManagerDefault & dof_manager;
};

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::updateDOFsData(
    DOFDataDefault & dof_data, UInt nb_new_local_dofs, UInt nb_new_pure_local,
    const std::function<UInt(UInt)> & getNode) {
  auto prank = this->communicator.whoAmI();
  auto psize = this->communicator.getNbProc();

  // access the relevant data to update
  const auto & support_type = dof_data.support_type;

  // resize all relevant arrays
  this->residual.resize(this->local_system_size, 0.);
  this->dofs_flag.resize(this->local_system_size, NodeFlag::_normal);
  this->global_solution.resize(this->local_system_size, 0.);
  this->global_blocked_dofs.resize(this->local_system_size, true);
  this->previous_global_blocked_dofs.resize(this->local_system_size, true);
  this->global_equation_number.resize(this->local_system_size, -1);

  for (auto & lumped_matrix : lumped_matrices)
    lumped_matrix.second->resize(this->local_system_size);

  dof_data.local_equation_number.reserve(dof_data.local_equation_number.size() +
                                         nb_new_local_dofs);

  // determine the first local/global dof id to use
  Array<UInt> nb_dofs_per_proc(psize);
  nb_dofs_per_proc(prank) = nb_new_pure_local;
  this->communicator.allGather(nb_dofs_per_proc);

  this->first_global_dof_id += std::accumulate(
      nb_dofs_per_proc.begin(), nb_dofs_per_proc.begin() + prank, 0);
  UInt first_dof_id = this->local_system_size - nb_new_local_dofs;

  if (support_type == _dst_nodal) {
    dof_data.associated_nodes.reserve(dof_data.associated_nodes.size() +
                                      nb_new_local_dofs);
  }

  // update per dof info
  for (UInt d = 0; d < nb_new_local_dofs; ++d) {
    UInt local_eq_num = first_dof_id + d;
    dof_data.local_equation_number.push_back(local_eq_num);

    bool is_local_dof = true;

    // determine the dof type
    switch (support_type) {
    case _dst_nodal: {
      UInt node = getNode(d / dof_data.dof->getNbComponent());

      this->dofs_flag(local_eq_num) = this->mesh->getNodeFlag(node);
      dof_data.associated_nodes.push_back(node);
      is_local_dof = this->mesh->isLocalOrMasterNode(node);
      break;
    }
    case _dst_generic: {
      this->dofs_flag(local_eq_num) = NodeFlag::_normal;
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
      this->global_equation_number(local_eq_num) = UInt(-1);
    }
  }

  // synchronize the global numbering for slaves
  if (support_type == _dst_nodal && this->synchronizer) {
    GlobalDOFInfoDataAccessor data_accessor(dof_data, *this);
    auto & node_synchronizer = this->mesh->getNodeSynchronizer();
    node_synchronizer.synchronizeOnce(data_accessor, _gst_ask_nodes);
  }

  this->first_global_dof_id += std::accumulate(
      nb_dofs_per_proc.begin() + prank + 1, nb_dofs_per_proc.end(), 0);

  matrix_profiled_dofs.clear();
  for(auto & matrix : matrices) {
    matrix.second->clearProfile();
  }

}

/* -------------------------------------------------------------------------- */
// register in factory
static bool default_dof_manager_is_registered[[gnu::unused]] =
    DefaultDOFManagerFactory::getInstance().registerAllocator(
        "default", [](const ID & id,
                      const MemoryID & mem_id) -> std::unique_ptr<DOFManager> {
          return std::make_unique<DOFManagerDefault>(id, mem_id);
        });

static bool dof_manager_is_registered[[gnu::unused]] =
    DOFManagerFactory::getInstance().registerAllocator(
        "default", [](Mesh & mesh, const ID & id,
                      const MemoryID & mem_id) -> std::unique_ptr<DOFManager> {
          return std::make_unique<DOFManagerDefault>(mesh, id, mem_id);
        });
} // namespace akantu
