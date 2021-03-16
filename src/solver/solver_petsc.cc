/**
 * @file   solver_petsc.cc
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue May 13 2014
 * @date last modification: Sun Aug 13 2017
 *
 * @brief  Solver class implementation for the petsc solver
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "solver_petsc.hh"
#include "dof_manager_petsc.hh"
#include "mpi_communicator_data.hh"
#include "solver_vector_petsc.hh"
#include "sparse_matrix_petsc.hh"
/* -------------------------------------------------------------------------- */
#include <petscksp.h>
//#include <petscsys.h>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
SolverPETSc::SolverPETSc(DOFManagerPETSc & dof_manager, const ID & matrix_id,
                         const ID & id)
    : SparseSolver(dof_manager, matrix_id, id),
      dof_manager(dof_manager), matrix(dof_manager.getMatrix(matrix_id)) {
  auto && mpi_comm = dof_manager.getMPIComm();

  /// create a solver context
  PETSc_call(KSPCreate, mpi_comm, &this->ksp);
}

/* -------------------------------------------------------------------------- */
SolverPETSc::~SolverPETSc() {
  if (ksp != nullptr) {
    PETSc_call(KSPDestroy, &ksp);
  }
}

/* -------------------------------------------------------------------------- */
void SolverPETSc::setOperators() {
  // set the matrix that defines the linear system and the matrix for
// preconditioning (here they are the same)
#if PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5
  PETSc_call(KSPSetOperators, ksp, this->matrix.getMat(),
             this->matrix.getMat());
#else
  PETSc_call(KSPSetOperators, ksp, this->matrix.getMat(), this->matrix.getMat(),
             SAME_NONZERO_PATTERN);
#endif

  // If this is not called the solution vector is zeroed in the call to
  // KSPSolve().
  PETSc_call(KSPSetInitialGuessNonzero, ksp, PETSC_TRUE);
  PETSc_call(KSPSetFromOptions, ksp);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolverPETSc::solve() {
  Vec & rhs(this->dof_manager.getResidual());
  Vec & solution(this->dof_manager.getSolution());

  PETSc_call(KSPSolve, ksp, rhs, solution);
}

} // namespace akantu
