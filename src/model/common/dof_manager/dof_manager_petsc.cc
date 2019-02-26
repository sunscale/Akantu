/**
 * @file   dof_manager_petsc.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Oct 07 2015
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  DOFManaterPETSc is the PETSc implementation of the DOFManager
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
#include "dof_manager_petsc.hh"
#include "communicator.hh"
#include "cppargparse.hh"
#include "non_linear_solver_petsc.hh"
#include "solver_vector_petsc.hh"
#include "sparse_matrix_petsc.hh"
#include "time_step_solver_default.hh"
#if defined(AKANTU_USE_MPI)
#include "mpi_communicator_data.hh"
#endif
/* -------------------------------------------------------------------------- */
#include <petscis.h>
#include <petscsys.h>
/* -------------------------------------------------------------------------- */

namespace akantu {

int DOFManagerPETSc::petsc_dof_manager_instances = 0;

/* -------------------------------------------------------------------------- */
namespace internal {
  PETScVector::~PETScVector() = default;

  PETScVector::operator Vec() { return x; }

  Int PETScVector::size() const {
    PetscInt n;
    PETSc_call(VecGetSize, x, &n);
    return n;
  }

  Int PETScVector::local_size() const {
    PetscInt n;
    PETSc_call(VecGetLocalSize, x, &n);
    return n;
  }

  template <class Array>
  PETScWrapedVector<Array>::PETScWrapedVector(Array && array) : array(array) {
    PETSc_call(VecCreateSeqWithArray, PETSC_COMM_SELF, 1, array.size(),
               array.storage(), &x);
  }

  template <class Array> PETScWrapedVector<Array>::~PETScWrapedVector() {
    PETSc_call(VecDestroy, &x);
  }

  template <class V>
  PETScLocalVector<V>::PETScLocalVector(const SolverVector & g)
      : g(dynamic_cast<const SolverVectorPETSc &>(g).getVec()) {
    PETSc_call(VecGetLocalVector, this->g, x);
  }

  template <class V> PETScLocalVector<V>::~PETScLocalVector() {
    PETSc_call(VecRestoreLocalVector, x, g);
    PETSc_call(VecDestroy, &x);
  }

  PETScLocalVector<const Vec &> make_petsc_local_vector(const SolverVector & vec) {
    return PETScLocalVector<const Vec &>(vec);
  }

} // namespace internal

/* -------------------------------------------------------------------------- */
DOFManagerPETSc::DOFDataPETSc::DOFDataPETSc(const ID & dof_id)
    : DOFData(dof_id) {}

/* -------------------------------------------------------------------------- */
DOFManagerPETSc::DOFManagerPETSc(const ID & id, const MemoryID & memory_id)
    : DOFManager(id, memory_id) {
  init();
}

/* -------------------------------------------------------------------------- */
DOFManagerPETSc::DOFManagerPETSc(Mesh & mesh, const ID & id,
                                 const MemoryID & memory_id)
    : DOFManager(mesh, id, memory_id) {
  init();
}

/* -------------------------------------------------------------------------- */
void DOFManagerPETSc::init() {
  // check if the akantu types and PETSc one are consistant
  static_assert(sizeof(Int) == sizeof(PetscInt),
                "The integer type of Akantu does not match the one from PETSc");
  static_assert(sizeof(Real) == sizeof(PetscReal),
                "The integer type of Akantu does not match the one from PETSc");

#if defined(AKANTU_USE_MPI)
  const auto & mpi_data =
      aka::as_type<MPICommunicatorData>(communicator.getCommunicatorData());
  MPI_Comm mpi_comm = mpi_data.getMPICommunicator();
  this->mpi_communicator = mpi_comm;
#else
  this->mpi_communicator = PETSC_COMM_SELF;
#endif

  PetscBool isInitialized;
  PETSc_call(PetscInitialized, &isInitialized);

  if (not isInitialized) {
    cppargparse::ArgumentParser & argparser = getStaticArgumentParser();
    int & argc = argparser.getArgC();
    char **& argv = argparser.getArgV();
    PETSc_call(PetscInitialize, &argc, &argv, NULL, NULL);
    PETSc_call(PetscPopErrorHandler); // remove the default PETSc signal handler
    PETSc_call(PetscPushErrorHandler, PetscIgnoreErrorHandler, NULL);
  }

  this->petsc_dof_manager_instances++;
}

/* -------------------------------------------------------------------------- */
DOFManagerPETSc::~DOFManagerPETSc() {
  this->petsc_dof_manager_instances--;

  for (auto && data : dofs) {
    auto & dof_data_petsc = dynamic_cast<DOFDataPETSc &>(*(data.second));
    if (dof_data_petsc.is) {
      PETSc_call(ISDestroy, &dof_data_petsc.is);
    }
  }

  if (is_ltog_map)
    PETSc_call(ISLocalToGlobalMappingDestroy, &is_ltog_map);

  if (this->petsc_dof_manager_instances == 0) {
    PetscFinalize();
  }
}

/* -------------------------------------------------------------------------- */
auto DOFManagerPETSc::getNewDOFData(const ID & dof_id)
    -> std::unique_ptr<DOFData> {
  return std::make_unique<DOFDataPETSc>(dof_id);
}

/* -------------------------------------------------------------------------- */
std::tuple<UInt, UInt, UInt>
DOFManagerPETSc::registerDOFsInternal(const ID & dof_id,
                                      Array<Real> & dofs_array) {
  dofs_ids.push_back(dof_id);
  auto ret = DOFManager::registerDOFsInternal(dof_id, dofs_array);
  UInt nb_dofs, nb_pure_local_dofs;
  std::tie(nb_dofs, nb_pure_local_dofs, std::ignore) = ret;

  auto & dof_data = this->getDOFDataTyped<DOFDataPETSc>(dof_id);

  ISCreateGeneral(PETSC_COMM_SELF, dof_data.local_equation_number.size(),
                  dof_data.local_equation_number.storage(), PETSC_USE_POINTER,
                  &(dof_data.is));

  if (is_ltog_map)
    PETSc_call(ISLocalToGlobalMappingDestroy, &is_ltog_map);

  PETSc_call(ISLocalToGlobalMappingCreate, PETSC_COMM_SELF, 1,
             global_equation_number.size(), global_equation_number.storage(),
             PETSC_COPY_VALUES, &is_ltog_map);

  return ret;
}

/* -------------------------------------------------------------------------- */
void DOFManagerPETSc::assembleToGlobalArray(
    const ID & dof_id, const Array<Real> & array_to_assemble,
    SolverVector & global_array, Real scale_factor) {
  const auto & is = getDOFDataTyped<DOFDataPETSc>(dof_id).is;

  auto y = internal::make_petsc_local_vector(global_array);
  auto x = internal::make_petsc_wraped_vector(array_to_assemble);

  PETSc_call(VecISAXPY, y, is, scale_factor, x);
}

/* -------------------------------------------------------------------------- */
void DOFManagerPETSc::getArrayPerDOFs(const ID & dof_id,
                                      const SolverVector & global,
                                      Array<Real> & local) {
  const auto & is = getDOFDataTyped<DOFDataPETSc>(dof_id).is;

  auto g = internal::make_petsc_local_vector(global);

  const PetscInt * idx;
  PETSc_call(ISGetIndices, is, &idx);

  PETSc_call(VecGetValues, g, local.size(), idx, local.storage());

  PETSc_call(ISRestoreIndices, is, &idx);
}

/* -------------------------------------------------------------------------- */
void DOFManagerPETSc::makeConsistentForPeriodicity(const ID & /*dof_id*/,
                                                   SolverVector & /*array*/) {

}



/* -------------------------------------------------------------------------- */
NonLinearSolver &
DOFManagerPETSc::getNewNonLinearSolver(const ID & id,
                                       const NonLinearSolverType & type) {
  return this->registerNonLinearSolver<NonLinearSolverPETSc>(*this, id, type);
}

/* -------------------------------------------------------------------------- */
TimeStepSolver &
DOFManagerPETSc::getNewTimeStepSolver(const ID & id,
                                      const TimeStepSolverType & type,
                                      NonLinearSolver & non_linear_solver) {
  return this->registerTimeStepSolver<TimeStepSolverDefault>(*this, id, type,
                                                             non_linear_solver);
}

/* -------------------------------------------------------------------------- */
SparseMatrix & DOFManagerPETSc::getNewMatrix(const ID & id,
                                             const MatrixType & matrix_type) {
  return this->registerSparseMatrix<SparseMatrixPETSc>(*this, id, matrix_type);
}

/* -------------------------------------------------------------------------- */
SparseMatrix & DOFManagerPETSc::getNewMatrix(const ID & id,
                                             const ID & matrix_to_copy_id) {
  return this->registerSparseMatrix<SparseMatrixPETSc>(id, matrix_to_copy_id);
}

/* -------------------------------------------------------------------------- */
SparseMatrixPETSc & DOFManagerPETSc::getMatrix(const ID & id) {
  auto & matrix = DOFManager::getMatrix(id);
  return dynamic_cast<SparseMatrixPETSc &>(matrix);
}

/* -------------------------------------------------------------------------- */
SolverVector & DOFManagerPETSc::getNewLumpedMatrix(const ID & id) {
  return this->registerLumpedMatrix<SolverVectorPETSc>(*this, id);
}

/* -------------------------------------------------------------------------- */
SolverVectorPETSc & DOFManagerPETSc::getSolution() {
  return dynamic_cast<SolverVectorPETSc &>(*this->solution);
}

const SolverVectorPETSc & DOFManagerPETSc::getSolution() const {
  return dynamic_cast<const SolverVectorPETSc &>(*this->solution);
}

SolverVectorPETSc & DOFManagerPETSc::getResidual() {
  return dynamic_cast<SolverVectorPETSc &>(*this->residual);
}

const SolverVectorPETSc & DOFManagerPETSc::getResidual() const {
  return dynamic_cast<const SolverVectorPETSc &>(*this->residual);
}

/* -------------------------------------------------------------------------- */
static bool dof_manager_is_registered [[gnu::unused]] =
    DOFManagerFactory::getInstance().registerAllocator(
        "petsc",
        [](Mesh & mesh, const ID & id,
           const MemoryID & mem_id) -> std::unique_ptr<DOFManager> {
          return std::make_unique<DOFManagerPETSc>(mesh, id, mem_id);
        });

} // namespace akantu
