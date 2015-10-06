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

#if defined(AKANTU_USE_MPI)
#  include "mpi_type_wrapper.hh"
#endif

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_PETSC
#include <petscsys.h>
#endif

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
StaticSolver::StaticSolver() : CommunicatorEventHandler(), is_initialized(false) {
  StaticCommunicator::getStaticCommunicator().registerEventHandler(*this);
}


/* -------------------------------------------------------------------------- */
StaticSolver::~StaticSolver() {
  --this->nb_references;
  if(this->nb_references == 0) {
    StaticCommunicator::getStaticCommunicator().unregisterEventHandler(*this);
    delete this->static_solver;
  }
}

/* -------------------------------------------------------------------------- */
StaticSolver & StaticSolver::getStaticSolver() {
  if(nb_references == 0)
    static_solver = new StaticSolver();
  ++nb_references;
  return *static_solver;
}

#ifdef AKANTU_USE_PETSC
#endif

/* -------------------------------------------------------------------------- */
void StaticSolver::initialize(int & argc, char ** & argv) {
  if (this->is_initialized) return;
  //  AKANTU_DEBUG_ASSERT(this->is_initialized != true, "The static solver has already been initialized");


  this->is_initialized = true;
}

/* -------------------------------------------------------------------------- */
void StaticSolver::finalize() {
  ParentEventHandler::sendEvent(StaticSolverEvent::BeforeStaticSolverDestroyEvent());


  AKANTU_DEBUG_ASSERT(this->is_initialized == true, "The static solver has not been initialized");
#ifdef AKANTU_USE_PETSC
  PetscFinalize();
#endif

  this->is_initialized = false;
}

__END_AKANTU__
