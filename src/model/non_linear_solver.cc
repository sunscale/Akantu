/**
 * @file   non_linear_solver.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Tue Oct 13 15:34:43 2015
 *
 * @brief  Implementation of the base class NonLinearSolver
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
#include "non_linear_solver.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
void NonLinearSolver::callbackPredictors() {
  NonLinearSolverCallbackMap it = solver_callbacks.begin();
  NonLinearSolverCallbackMap end = solver_callbacks.end();

  for (; it != end; ++it) {
    it->second->predictor();
  }
}

/* -------------------------------------------------------------------------- */
void NonLinearSolver::callbackCorrectors() {
  NonLinearSolverCallbackMap it = solver_callbacks.begin();
  NonLinearSolverCallbackMap end = solver_callbacks.end();

  for (; it != end; ++it) {
    it->second->corrector();
  }
}

/* -------------------------------------------------------------------------- */
void NonLinearSolver::callbackAssembleJacobian() {
  NonLinearSolverCallbackMap it = solver_callbacks.begin();
  NonLinearSolverCallbackMap end = solver_callbacks.end();

  for (; it != end; ++it) {
    it->second->assembleJacobian();
  }
}

/* -------------------------------------------------------------------------- */
void NonLinearSolver::callbackAssembleResidual() {
  NonLinearSolverCallbackMap it = solver_callbacks.begin();
  NonLinearSolverCallbackMap end = solver_callbacks.end();

  for (; it != end; ++it) {
    it->second->assembleResidual();
  }
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__
