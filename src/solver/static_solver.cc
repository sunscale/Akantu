/**
 * @file   static_solver.cc
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Wed Jul 30 15:35:01 2014
 *
 * @brief  implementation of the static solver
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
#include "static_solver.hh"
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_PETSC
#include <petscsys.h>
#endif

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
StaticSolver::~StaticSolver() {
  --this->nb_references;
  if(this->nb_references == 0)
    delete this->static_solver;
}

/* -------------------------------------------------------------------------- */
StaticSolver & StaticSolver::getStaticSolver() {
  if(nb_references == 0)
    static_solver = new StaticSolver();
  
  ++nb_references;
  
  return *static_solver;
}

#ifdef AKANTU_USE_PETSC
#if PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5
static PetscErrorCode PETScErrorHandler(MPI_Comm,
					int line, const char * dir, const char *file,
					PetscErrorCode number,
					PetscErrorType type,
					const char *message,
					void *) {
  AKANTU_DEBUG_ERROR("An error occured in PETSc in file \"" << file << ":" << line << "\" - PetscErrorCode "<< number << " - \""<< message << "\"");
}
#else
static PetscErrorCode PETScErrorHandler(MPI_Comm,
					int line, const char * func, const char * dir, const char *file,
					PetscErrorCode number,
					PetscErrorType type,
					const char *message,
					void *) {
  AKANTU_DEBUG_ERROR("An error occured in PETSc in file \"" << file << ":" << line << "\" - PetscErrorCode "<< number << " - \""<< message << "\"");
}
#endif
#endif

/* -------------------------------------------------------------------------- */
void StaticSolver::initialize(int & argc, char ** & argv) {
  AKANTU_DEBUG_ASSERT(this->is_initialized != true, "The static solver has already been initialized");
#ifdef AKANTU_USE_PETSC
    PetscErrorCode petsc_error = PetscInitialize(&argc, &argv, NULL, NULL);
    if(petsc_error != 0) {
      AKANTU_DEBUG_ERROR("An error occured while initializing Petsc (PetscErrorCode "<< petsc_error << ")");
    }
    PetscPushErrorHandler(PETScErrorHandler, NULL);
#endif

    this->is_initialized = true;
  }

/* -------------------------------------------------------------------------- */
void StaticSolver::finalize() {
  AKANTU_DEBUG_ASSERT(this->is_initialized == true, "The static solver has not been initialized");
#ifdef AKANTU_USE_PETSC
  PetscFinalize();
#endif

  this->is_initialized = false;
}

__END_AKANTU__
