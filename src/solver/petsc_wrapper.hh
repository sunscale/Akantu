/**
 * @file   petsc_wrapper.hh
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Mon Jul 21 17:40:41 2014
 *
 * @brief  Wrapper of PETSc structures
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

#ifndef __AKANTU_PETSC_WRAPPER_HH__
#define __AKANTU_PETSC_WRAPPER_HH__

/* -------------------------------------------------------------------------- */
#include <mpi.h>
#include <petscmat.h>
#include <petscao.h>
#include <petscis.h>
#include <petscksp.h>

__BEGIN_AKANTU__

struct PETScMatrixWrapper {
  Mat mat;
  AO ao;
  ISLocalToGlobalMapping mapping;
  /// pointer to the MPI communicator for PETSc commands
  MPI_Comm communicator;
};

struct PETScSolverWrapper {
  KSP ksp;
  Vec solution;
  Vec rhs;
  /// pointer to the MPI communicator for PETSc commands
  MPI_Comm communicator;
};


__END_AKANTU__

#endif /* __AKANTU_PETSC_WRAPPER_HH__ */
