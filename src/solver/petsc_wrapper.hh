/**

 * @file   petsc_wrapper.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Feb 21 2013
 * @date last modification: Sat Feb 03 2018
 *
 * @brief  Wrapper of PETSc structures
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

#ifndef __AKANTU_PETSC_WRAPPER_HH__
#define __AKANTU_PETSC_WRAPPER_HH__

/* -------------------------------------------------------------------------- */
#include <mpi.h>
#include <petscao.h>
#include <petscis.h>
#include <petscksp.h>
#include <petscmat.h>
#include <petscvec.h>

namespace akantu {

/* -------------------------------------------------------------------------- */
struct PETScMatrixWrapper {
  Mat mat;
  AO ao;
  ISLocalToGlobalMapping mapping;
  /// MPI communicator for PETSc commands
  MPI_Comm communicator;
};

/* -------------------------------------------------------------------------- */
struct PETScSolverWrapper {
  KSP ksp;
  Vec solution;
  Vec rhs;
  // MPI communicator for PETSc commands
  MPI_Comm communicator;
};

#if not defined(PETSC_CLANGUAGE_CXX)
extern int aka_PETScError(int ierr);

#define CHKERRXX(x)                                                            \
  do {                                                                         \
    int error = aka_PETScError(x);                                             \
    if (error != 0) {                                                          \
      AKANTU_EXCEPTION("Error in PETSC");                                      \
    }                                                                          \
  } while (0)
#endif

} // namespace akantu

#endif /* __AKANTU_PETSC_WRAPPER_HH__ */
