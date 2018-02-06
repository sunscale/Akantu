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
#include "sparse_matrix_petsc.hh"
#include "non_linear_solver_petsc.hh"
#include "time_step_solver_default.hh"
#if defined(AKANTU_USE_MPI)
#include "mpi_communicator_data.hh"
#endif
/* -------------------------------------------------------------------------- */
#include <petscsys.h>
/* -------------------------------------------------------------------------- */

namespace akantu {

int DOFManagerPETSc::petsc_dof_manager_instances = 0;

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
#if __cplusplus > 199711L
  static_assert(sizeof(Int) == sizeof(PetscInt),
                "The integer type of Akantu does not match the one from PETSc");
  static_assert(sizeof(Real) == sizeof(PetscReal),
                "The integer type of Akantu does not match the one from PETSc");
#else
  AKANTU_DEBUG_ASSERT(
      sizeof(Int) == sizeof(PetscInt),
      "The integer type of Akantu does not match the one from PETSc");
  AKANTU_DEBUG_ASSERT(
      sizeof(Real) == sizeof(PetscReal),
      "The integer type of Akantu does not match the one from PETSc");
#endif

#if defined(AKANTU_USE_MPI)
  const auto & mpi_data = dynamic_cast<const MPICommunicatorData &>(
      communicator.getCommunicatorData());
  MPI_Comm mpi_comm = mpi_data.getMPICommunicator();
  this->mpi_communicator = mpi_comm;
#else
  this->mpi_communicator = PETSC_COMM_SELF;
#endif

  PetscBool isInitialized;
  PETSc_call(PetscInitialized, &isInitialized);
  if (isInitialized) {
    cppargparse::ArgumentParser & argparser = getStaticArgumentParser();
    int & argc = argparser.getArgC();
    char **& argv = argparser.getArgV();
    PETSc_call(PetscInitialize, &argc, &argv, NULL, NULL);
    PETSc_call(PetscPushErrorHandler, PetscIgnoreErrorHandler, NULL);
  }

  this->petsc_dof_manager_instances++;

  PETSc_call(VecCreate, mpi_communicator, &this->residual);
  PETSc_call(VecCreate, mpi_communicator, &this->solution);
}

/* -------------------------------------------------------------------------- */
DOFManagerPETSc::~DOFManagerPETSc() {
  PETSc_call(VecDestroy, &(this->residual));
  PETSc_call(VecDestroy, &(this->solution));

  this->petsc_dof_manager_instances--;
  if (this->petsc_dof_manager_instances == 0) {
    PetscFinalize();
  }
}

/* -------------------------------------------------------------------------- */
void DOFManagerPETSc::registerDOFs(const ID & dof_id, Array<Real> & dofs_array,
                                   DOFSupportType & support_type) {
  DOFManager::registerDOFs(dof_id, dofs_array, support_type);

  PetscInt current_size;
  PETSc_call(VecGetSize, this->residual, &current_size);

  PetscInt local_size = this->pure_local_system_size;
  PETSc_call(VecSetSizes, this->residual, local_size, system_size);

  PetscInt start_dof, end_dof;
  VecGetOwnershipRange(this->residual, &start_dof, &end_dof);

  std::vector<PetscInt> global_indices(local_size);
  global_indices[0] = start_dof;

  for (PetscInt d = 0; d < local_size; d++)
    global_indices[d + 1] = global_indices[d] + 1;

// To be change if we switch to a block definition
#if PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5
  ISLocalToGlobalMappingCreate(this->mpi_communicator, 1, local_size,
                               global_indices.data(), PETSC_COPY_VALUES,
                               &this->is_ltog_map);

#else
  ISLocalToGlobalMappingCreate(this->mpi_communicator, local_size,
                               global_indices, PETSC_COPY_VALUES,
                               &this->is_ltog_map);
#endif

  VecSetLocalToGlobalMapping(this->residual, this->is_ltog_map);
}

/* -------------------------------------------------------------------------- */
void DOFManagerPETSc::assembleToResidual(const ID & dof_id,
                                         const Array<Real> & array_to_assemble,
                                         Real scale_factor) {
  const auto & is = getDOFDataTyped<DOFDataPETSc>(dof_id).is;
  Vec y;
  PETSc_call(VecGetLocalVector, this->residual, y);

  Vec x;
  PETSc_call(VecCreateSeqWithArray, PETSC_COMM_SELF, 1,
             array_to_assemble.size(), array_to_assemble.storage(), &x);
  PETSc_call(VecISAXPY, y, is, scale_factor, x);
  PETSc_call(VecRestoreLocalVector, y, this->residual);
  PETSc_call(VecDestroy, &y);
  PETSc_call(VecDestroy, &x);
}

/* -------------------------------------------------------------------------- */
void DOFManagerPETSc::clearResidual() {
  PETSc_call(VecZeroEntries, this->residual);
}

/* -------------------------------------------------------------------------- */
void DOFManagerPETSc::clearMatrix(const ID & mtx) {
  auto & matrix = this->getMatrix(mtx);
  matrix.clear();
}

/* -------------------------------------------------------------------------- */
void DOFManagerPETSc::clearLumpedMatrix(const ID & mtx) {
  auto & global_lumped = this->getLumpedMatrix(mtx);
  PETSc_call(VecZeroEntries, global_lumped);
}

/* -------------------------------------------------------------------------- */
static void getLocalValues(Vec & garray, Array<Real> & larray,
                           Array<Int> & lidx) {
  PetscInt n = lidx.size();
  Array<PetscInt> global_idx(n);

  ISLocalToGlobalMapping mapping;
  PETSc_call(VecGetLocalToGlobalMapping, garray, &mapping);
  PETSc_call(ISLocalToGlobalMappingApply, mapping, n, lidx.storage(),
             global_idx.storage());

  PETSc_call(VecAssemblyBegin, garray);
  PETSc_call(VecAssemblyEnd, garray);
  PETSc_call(VecGetValues, garray, n, global_idx.storage(), larray.storage());
}

/* -------------------------------------------------------------------------- */
void DOFManagerPETSc::getLumpedMatrixPerDOFs(const ID & dof_id,
                                             const ID & lumped_mtx,
                                             Array<Real> & lumped) {
  auto & dof_data = this->getDOFData(dof_id);
  auto & global_lumped = this->getLumpedMatrix(lumped_mtx);
  getLocalValues(global_lumped, lumped, dof_data.local_equation_number);
}

/* -------------------------------------------------------------------------- */
void DOFManagerPETSc::getSolutionPerDOFs(const ID & dof_id,
                                         Array<Real> & solution_array) {
  auto & dof_data = this->getDOFData(dof_id);
  getLocalValues(this->solution, solution_array,
                 dof_data.local_equation_number);
}

/* -------------------------------------------------------------------------- */
NonLinearSolver & DOFManagerPETSc::getNewNonLinearSolver(
    const ID & id,
    const NonLinearSolverType & type) {
  ID non_linear_solver_id = this->id + ":nls:" + id;
  std::unique_ptr<NonLinearSolver> nls = std::make_unique<NonLinearSolverPETSc>(
        *this, type, non_linear_solver_id, this->memory_id);

  return this->registerNonLinearSolver(non_linear_solver_id, nls);
}

/* -------------------------------------------------------------------------- */
TimeStepSolver &
DOFManagerPETSc::getNewTimeStepSolver(const ID & id,
                                      const TimeStepSolverType & type,
                                      NonLinearSolver & non_linear_solver) {
  ID time_step_solver_id = this->id + ":tss:" + id;

  std::unique_ptr<TimeStepSolver> tss = std::make_unique<TimeStepSolverDefault>(
      *this, type, non_linear_solver, time_step_solver_id, this->memory_id);

  return this->registerTimeStepSolver(time_step_solver_id, tss);
}

/* -------------------------------------------------------------------------- */
static bool dof_manager_is_registered[[gnu::unused]] =
    DOFManagerFactory::getInstance().registerAllocator(
        "petsc", [](Mesh & mesh, const ID & id,
                    const MemoryID & mem_id) -> std::unique_ptr<DOFManager> {
          return std::make_unique<DOFManagerPETSc>(mesh, id, mem_id);
        });

} // akantu
