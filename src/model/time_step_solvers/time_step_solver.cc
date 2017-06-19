/**
 * @file   time_step_solver.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Mon Oct 12 16:56:43 2015
 *
 * @brief  Implementation of common part of TimeStepSolvers
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "time_step_solver.hh"
#include "non_linear_solver.hh"
#include "dof_manager.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
TimeStepSolver::TimeStepSolver(DOFManager & dof_manager,
                               const TimeStepSolverType & type,
                               NonLinearSolver & non_linear_solver,
                               const ID & id, UInt memory_id)
    : Memory(id, memory_id), SolverCallback(dof_manager),
      _dof_manager(dof_manager), type(type), time_step(0.),
      solver_callback(NULL), non_linear_solver(non_linear_solver) {
  this->registerSubRegistry("non_linear_solver", non_linear_solver);
}

/* -------------------------------------------------------------------------- */
TimeStepSolver::~TimeStepSolver() {}

/* -------------------------------------------------------------------------- */
void TimeStepSolver::predictor() {
  AKANTU_DEBUG_ASSERT(
      this->solver_callback != NULL,
      "This function cannot be called if the solver_callback is not set");

  this->solver_callback->predictor();
}

/* -------------------------------------------------------------------------- */
void TimeStepSolver::corrector() {
  AKANTU_DEBUG_ASSERT(
      this->solver_callback != NULL,
      "This function cannot be called if the solver_callback is not set");

  this->solver_callback->corrector();
}

/* -------------------------------------------------------------------------- */
void TimeStepSolver::assembleJacobian() {
  AKANTU_DEBUG_ASSERT(
      this->solver_callback != NULL,
      "This function cannot be called if the solver_callback is not set");

  //  this->_dof_manager.clearMatrix("J");
  this->solver_callback->assembleJacobian();
}

/* -------------------------------------------------------------------------- */
void TimeStepSolver::assembleResidual() {
  AKANTU_DEBUG_ASSERT(
      this->solver_callback != NULL,
      "This function cannot be called if the solver_callback is not set");

  this->_dof_manager.clearResidual();
  this->solver_callback->assembleResidual();
}

/* -------------------------------------------------------------------------- */

} // akantu
