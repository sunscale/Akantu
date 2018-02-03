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

#include "cppargparse.hh"

#if defined(AKANTU_USE_MPI)
#include "mpi_type_wrapper.hh"
#include "static_communicator.hh"
#endif
/* -------------------------------------------------------------------------- */
#include <petscsys.h>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
DOFManagerPETSc::DOFManagerPETSc(const ID & id, const MemoryID & memory_id)
    : DOFManager(id, memory_id) {
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
  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  const StaticCommunicatorMPI & mpi_st_comm =
      dynamic_cast<const StaticCommunicatorMPI &>(
          comm.getRealStaticCommunicator());

  this->mpi_communicator = mpi_st_comm.getMPITypeWrapper().getMPICommunicator();
#else
  this->mpi_communicator = PETSC_COMM_SELF;
#endif

  PetscBool isInitialized;
  PETSc_call(PetscInitialized, &isInitialized);
  if(isInitialized) {
    cppargparse::ArgumentParser & argparser = getStaticArgumentParser();
    int & argc = argparser.getArgC();
    char **& argv = argparser.getArgV();
    PETSc_call(PetscInitialize, &argc, &argv, NULL, NULL);
    PETSc_call(PetscPushErrorHandler, PETScErrorHandler, NULL);
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

  if (current_size == 0) { // first time vector is set
    PetscInt local_size = this->pure_local_system_size;
    PETSc_call(VecSetSizes, this->residual, local_size, PETSC_DECIDE);
    PETSc_call(VecSetFromOptions, this->residual);

#ifndef AKANTU_NDEBUG
    PetscInt global_size;
    ierr = VecGetSize(this->residual, &global_size);
    CHKERRXX(ierr);
    AKANTU_DEBUG_ASSERT(this->system_size == UInt(global_size),
                        "The local value of the system size does not match the "
                        "one determined by PETSc");
#endif
    PetscInt start_dof, end_dof;
    VecGetOwnershipRange(this->residual, &start_dof, &end_dof);

    PetscInt * global_indices = new PetscInt[local_size];
    global_indices[0] = start_dof;

    for (PetscInt d = 0; d < local_size; d++)
      global_indices[d + 1] = global_indices[d] + 1;

// To be change if we switch to a block definition
#if PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5
    ISLocalToGlobalMappingCreate(this->communicator, 1, local_size,
                                 global_indices, PETSC_COPY_VALUES,
                                 &this->is_ltog);

#else
    ISLocalToGlobalMappingCreate(this->communicator, local_size, global_indices,
                                 PETSC_COPY_VALUES, &this->is_ltog);
#endif

    VecSetLocalToGlobalMapping(this->residual, this->is_ltog);
    delete[] global_indices;

    ierr = VecDuplicate(this->residual, &this->solution);
    CHKERRXX(ierr);

  } else { // this is an update of the object already created
    AKANTU_TO_IMPLEMENT();
  }

  /// set the solution to zero
  // ierr = VecZeroEntries(this->solution);
  // CHKERRXX(ierr);
}

} // akantu
