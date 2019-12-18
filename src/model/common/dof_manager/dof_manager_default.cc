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
#include "non_linear_solver_default.hh"
#include "periodic_node_synchronizer.hh"
#include "solver_vector_default.hh"
#include "solver_vector_distributed.hh"
#include "sparse_matrix_aij.hh"
#include "time_step_solver_default.hh"

/* -------------------------------------------------------------------------- */
#include <algorithm>
#include <memory>
#include <numeric>
#include <unordered_map>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
DOFManagerDefault::DOFManagerDefault(const ID & id, const MemoryID & memory_id)
    : DOFManager(id, memory_id), synchronizer(nullptr) {
  residual = std::make_unique<SolverVectorDefault>(
      *this, std::string(id + ":residual"));
  solution = std::make_unique<SolverVectorDefault>(
      *this, std::string(id + ":solution"));
  data_cache = std::make_unique<SolverVectorDefault>(
      *this, std::string(id + ":data_cache"));
}

/* -------------------------------------------------------------------------- */
DOFManagerDefault::DOFManagerDefault(Mesh & mesh, const ID & id,
                                     const MemoryID & memory_id)
    : DOFManager(mesh, id, memory_id), synchronizer(nullptr) {
  if (this->mesh->isDistributed()) {
    this->synchronizer = std::make_unique<DOFSynchronizer>(
        *this, this->id + ":dof_synchronizer", this->memory_id);
    residual = std::make_unique<SolverVectorDistributed>(
        *this, std::string(id + ":residual"));
    solution = std::make_unique<SolverVectorDistributed>(
        *this, std::string(id + ":solution"));
    data_cache = std::make_unique<SolverVectorDistributed>(
        *this, std::string(id + ":data_cache"));

  } else {
    residual = std::make_unique<SolverVectorDefault>(
        *this, std::string(id + ":residual"));
    solution = std::make_unique<SolverVectorDefault>(
        *this, std::string(id + ":solution"));
    data_cache = std::make_unique<SolverVectorDefault>(
        *this, std::string(id + ":data_cache"));
  }
}

/* -------------------------------------------------------------------------- */
DOFManagerDefault::~DOFManagerDefault() = default;

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::makeConsistentForPeriodicity(const ID & dof_id,
                                                     SolverVector & array) {
  auto & dof_data = this->getDOFDataTyped<DOFDataDefault>(dof_id);
  if (dof_data.support_type != _dst_nodal)
    return;

  if (not mesh->isPeriodic())
    return;

  this->mesh->getPeriodicNodeSynchronizer()
      .reduceSynchronizeWithPBCSlaves<AddOperation>(
          aka::as_type<SolverVectorDefault>(array).getVector());
}

/* -------------------------------------------------------------------------- */
template <typename T>
void DOFManagerDefault::assembleToGlobalArray(
    const ID & dof_id, const Array<T> & array_to_assemble,
    Array<T> & global_array, T scale_factor) {
  AKANTU_DEBUG_IN();

  auto & dof_data = this->getDOFDataTyped<DOFDataDefault>(dof_id);
  AKANTU_DEBUG_ASSERT(dof_data.local_equation_number.size() ==
                          array_to_assemble.size() *
                              array_to_assemble.getNbComponent(),
                      "The array to assemble does not have a correct size."
                          << " (" << array_to_assemble.getID() << ")");

  if (dof_data.support_type == _dst_nodal and mesh->isPeriodic()) {
    for (auto && data :
         zip(dof_data.local_equation_number, dof_data.associated_nodes,
             make_view(array_to_assemble))) {
      auto && equ_num = std::get<0>(data);
      auto && node = std::get<1>(data);
      auto && arr = std::get<2>(data);
      global_array(equ_num) +=
          scale_factor * (arr) * (not this->mesh->isPeriodicSlave(node));
    }
  } else {
    for (auto && data :
         zip(dof_data.local_equation_number, make_view(array_to_assemble))) {
      auto && equ_num = std::get<0>(data);
      auto && arr = std::get<1>(data);
      global_array(equ_num) += scale_factor * (arr);
    }
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::assembleToGlobalArray(
    const ID & dof_id, const Array<Real> & array_to_assemble,
    SolverVector & global_array_v, Real scale_factor) {

  assembleToGlobalArray(
      dof_id, array_to_assemble,
      aka::as_type<SolverVectorDefault>(global_array_v).getVector(),
      scale_factor);
}

/* -------------------------------------------------------------------------- */
DOFManagerDefault::DOFDataDefault::DOFDataDefault(const ID & dof_id)
    : DOFData(dof_id) {}

/* -------------------------------------------------------------------------- */
auto DOFManagerDefault::getNewDOFData(const ID & dof_id)
    -> std::unique_ptr<DOFData> {
  return std::make_unique<DOFDataDefault>(dof_id);
}

/* -------------------------------------------------------------------------- */
std::tuple<UInt, UInt, UInt>
DOFManagerDefault::registerDOFsInternal(const ID & dof_id,
                                        Array<Real> & dofs_array) {
  auto ret = DOFManager::registerDOFsInternal(dof_id, dofs_array);

  // update the synchronizer if needed
  if (this->synchronizer)
    this->synchronizer->registerDOFs(dof_id);

  return ret;
}

/* -------------------------------------------------------------------------- */
SparseMatrix & DOFManagerDefault::getNewMatrix(const ID & id,
                                               const MatrixType & matrix_type) {
  return this->registerSparseMatrix<SparseMatrixAIJ>(*this, id, matrix_type);
}

/* -------------------------------------------------------------------------- */
SparseMatrix & DOFManagerDefault::getNewMatrix(const ID & id,
                                               const ID & matrix_to_copy_id) {
  return this->registerSparseMatrix<SparseMatrixAIJ>(id, matrix_to_copy_id);
}

/* -------------------------------------------------------------------------- */
SolverVector & DOFManagerDefault::getNewLumpedMatrix(const ID & id) {
  return this->registerLumpedMatrix<SolverVectorDefault>(*this, id);
}

/* -------------------------------------------------------------------------- */
SparseMatrixAIJ & DOFManagerDefault::getMatrix(const ID & id) {
  auto & matrix = DOFManager::getMatrix(id);
  return aka::as_type<SparseMatrixAIJ>(matrix);
}

/* -------------------------------------------------------------------------- */
NonLinearSolver &
DOFManagerDefault::getNewNonLinearSolver(const ID & id,
                                         const NonLinearSolverType & type) {
  switch (type) {
#if defined(AKANTU_USE_MUMPS)
  case NonLinearSolverType::_newton_raphson:
    /* FALLTHRU */
    /* [[fallthrough]]; un-comment when compiler will get it */
  case NonLinearSolverType::_newton_raphson_modified: {
    return this->registerNonLinearSolver<NonLinearSolverNewtonRaphson>(
        *this, id, type);
  }
  case NonLinearSolverType::_linear: {
    return this->registerNonLinearSolver<NonLinearSolverLinear>(*this, id,
                                                                type);
  }
#endif
  case NonLinearSolverType::_lumped: {
    return this->registerNonLinearSolver<NonLinearSolverLumped>(*this, id,
                                                                type);
  }
  default:
    AKANTU_EXCEPTION("The asked type of non linear solver is not supported by "
                     "this dof manager");
  }
}

/* -------------------------------------------------------------------------- */
TimeStepSolver & DOFManagerDefault::getNewTimeStepSolver(
    const ID & id, const TimeStepSolverType & type,
    NonLinearSolver & non_linear_solver, SolverCallback & solver_callback) {
  return this->registerTimeStepSolver<TimeStepSolverDefault>(
      *this, id, type, non_linear_solver, solver_callback);
}

/* -------------------------------------------------------------------------- */
template <typename T>
void DOFManagerDefault::getArrayPerDOFs(const ID & dof_id,
                                        const Array<T> & global_array,
                                        Array<T> & local_array) const {
  AKANTU_DEBUG_IN();

  const Array<Int> & equation_number = this->getLocalEquationsNumbers(dof_id);

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
void DOFManagerDefault::getArrayPerDOFs(const ID & dof_id,
                                        const SolverVector & global_array,
                                        Array<Real> & local_array) {
  getArrayPerDOFs(dof_id,
                  aka::as_type<SolverVectorDefault>(global_array).getVector(),
                  local_array);
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::assembleLumpedMatMulVectToResidual(
    const ID & dof_id, const ID & A_id, const Array<Real> & x,
    Real scale_factor) {
  const Array<Real> & A = this->getLumpedMatrix(A_id);
  auto & cache = aka::as_type<SolverVectorArray>(*this->data_cache);

  cache.clear();
  this->assembleToGlobalArray(dof_id, x, cache.getVector(), scale_factor);

  for (auto && data : zip(make_view(A), make_view(cache.getVector()),
                          make_view(this->getResidualArray()))) {
    const auto & A = std::get<0>(data);
    const auto & x = std::get<1>(data);
    auto & r = std::get<2>(data);
    r += A * x;
  }
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::assembleElementalMatricesToMatrix(
    const ID & matrix_id, const ID & dof_id, const Array<Real> & elementary_mat,
    const ElementType & type, const GhostType & ghost_type,
    const MatrixType & elemental_matrix_type,
    const Array<UInt> & filter_elements) {
  this->addToProfile(matrix_id, dof_id, type, ghost_type);
  auto & A = getMatrix(matrix_id);
  DOFManager::assembleElementalMatricesToMatrix_(
      A, dof_id, elementary_mat, type, ghost_type, elemental_matrix_type,
      filter_elements);
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::assemblePreassembledMatrix(
    const ID & dof_id_m, const ID & dof_id_n, const ID & matrix_id,
    const TermsToAssemble & terms) {
  auto & A = getMatrix(matrix_id);
  DOFManager::assemblePreassembledMatrix_(A, dof_id_m, dof_id_n, terms);
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::assembleMatMulVectToArray(const ID & dof_id,
                                                  const ID & A_id,
                                                  const Array<Real> & x,
                                                  Array<Real> & array,
                                                  Real scale_factor) {
  if (mesh->isDistributed()) {
    DOFManager::assembleMatMulVectToArray_<SolverVectorDistributed>(
        dof_id, A_id, x, array, scale_factor);
  } else {
    DOFManager::assembleMatMulVectToArray_<SolverVectorDefault>(
        dof_id, A_id, x, array, scale_factor);
  }
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::addToProfile(const ID & matrix_id, const ID & dof_id,
                                     const ElementType & type,
                                     const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  const auto & dof_data = this->getDOFData(dof_id);

  if (dof_data.support_type != _dst_nodal)
    return;

  auto mat_dof = std::make_pair(matrix_id, dof_id);
  auto type_pair = std::make_pair(type, ghost_type);

  auto prof_it = this->matrix_profiled_dofs.find(mat_dof);
  if (prof_it != this->matrix_profiled_dofs.end() &&
      std::find(prof_it->second.begin(), prof_it->second.end(), type_pair) !=
          prof_it->second.end())
    return;

  auto nb_degree_of_freedom_per_node = dof_data.dof->getNbComponent();

  const auto & equation_number = this->getLocalEquationsNumbers(dof_id);

  auto & A = this->getMatrix(matrix_id);
  A.resize(system_size);
  auto size = A.size();

  auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

  const auto & connectivity = this->mesh->getConnectivity(type, ghost_type);
  auto cbegin = connectivity.begin(nb_nodes_per_element);
  auto cit = cbegin;

  auto nb_elements = connectivity.size();
  UInt * ge_it = nullptr;
  if (dof_data.group_support != "__mesh__") {
    const auto & group_elements =
        this->mesh->getElementGroup(dof_data.group_support)
            .getElements(type, ghost_type);
    ge_it = group_elements.storage();
    nb_elements = group_elements.size();
  }

  UInt size_mat = nb_nodes_per_element * nb_degree_of_freedom_per_node;
  Vector<Int> element_eq_nb(size_mat);

  for (UInt e = 0; e < nb_elements; ++e) {
    if (ge_it)
      cit = cbegin + *ge_it;

    this->extractElementEquationNumber(
        equation_number, *cit, nb_degree_of_freedom_per_node, element_eq_nb);
    std::transform(
        element_eq_nb.storage(), element_eq_nb.storage() + element_eq_nb.size(),
        element_eq_nb.storage(),
        [&](auto & local) { return this->localToGlobalEquationNumber(local); });

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
Array<Real> & DOFManagerDefault::getSolutionArray() {
  return dynamic_cast<SolverVectorDefault *>(this->solution.get())->getVector();
}

/* -------------------------------------------------------------------------- */
const Array<Real> & DOFManagerDefault::getResidualArray() const {
  return dynamic_cast<SolverVectorDefault *>(this->residual.get())->getVector();
}

/* -------------------------------------------------------------------------- */
Array<Real> & DOFManagerDefault::getResidualArray() {
  return dynamic_cast<SolverVectorDefault *>(this->residual.get())->getVector();
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::onNodesAdded(const Array<UInt> & nodes_list,
                                     const NewNodesEvent & event) {
  DOFManager::onNodesAdded(nodes_list, event);

  if (this->synchronizer)
    this->synchronizer->onNodesAdded(nodes_list);
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::resizeGlobalArrays() {
  DOFManager::resizeGlobalArrays();

  this->global_blocked_dofs.resize(this->local_system_size, true);
  this->previous_global_blocked_dofs.resize(this->local_system_size, true);

  matrix_profiled_dofs.clear();
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::updateGlobalBlockedDofs() {
  DOFManager::updateGlobalBlockedDofs();

  if (this->global_blocked_dofs_release ==
      this->previous_global_blocked_dofs_release)
    return;

  global_blocked_dofs_uint.resize(local_system_size);
  global_blocked_dofs_uint.set(false);
  for (const auto & dof : global_blocked_dofs) {
    global_blocked_dofs_uint[dof] = true;
  }
}

/* -------------------------------------------------------------------------- */
Array<bool> & DOFManagerDefault::getBlockedDOFs() {
  return global_blocked_dofs_uint;
}

/* -------------------------------------------------------------------------- */
const Array<bool> & DOFManagerDefault::getBlockedDOFs() const {
  return global_blocked_dofs_uint;
}

/* -------------------------------------------------------------------------- */
// register in factory
// static bool default_dof_manager_is_registered [[gnu::unused]] =
//     DefaultDOFManagerFactory::getInstance().registerAllocator(
//         "default",
//         [](const ID & id,
//            const MemoryID & mem_id) -> std::unique_ptr<DOFManager> {
//           return std::make_unique<DOFManagerDefault>(id, mem_id);
//         });

static bool dof_manager_is_registered [[gnu::unused]] =
    DOFManagerFactory::getInstance().registerAllocator(
        "default",
        [](Mesh & mesh, const ID & id,
           const MemoryID & mem_id) -> std::unique_ptr<DOFManager> {
          return std::make_unique<DOFManagerDefault>(mesh, id, mem_id);
        });

static bool dof_manager_is_registered_mumps [[gnu::unused]] =
    DOFManagerFactory::getInstance().registerAllocator(
        "mumps",
        [](Mesh & mesh, const ID & id,
           const MemoryID & mem_id) -> std::unique_ptr<DOFManager> {
          return std::make_unique<DOFManagerDefault>(mesh, id, mem_id);
        });

} // namespace akantu
