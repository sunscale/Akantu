/**
 * @file   time_step_solver_default_solver_callback.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Mon Sep 28 18:58:07 2015
 *
 * @brief  Implementation of the NonLinearSolverCallback for the time step
 * solver
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
#include "non_linear_solver_callback.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_TIME_STEP_SOLVER_DEFAULT_SOLVER_CALLBACK_HH__
#define __AKANTU_TIME_STEP_SOLVER_DEFAULT_SOLVER_CALLBACK_HH__

namespace akantu {

class TimeStepSolverCallback : public NonLinearSolverCallback {
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// callback to assemble the Jacobian Matrix
  virtual void assembleJacobian();

  /// callback to assemble the residual (rhs)
  virtual void assembleResidual();

  /* ------------------------------------------------------------------------ */
  /* Dynamic simulations part                                                 */
  /* ------------------------------------------------------------------------ */
  /// callback for the predictor (in case of dynamic simulation)
  virtual void predictor();

  /// callback for the corrector (in case of dynamic simulation)
  virtual void corrector();
};

} // akantu

#endif /* __AKANTU_TIME_STEP_SOLVER_DEFAULT_SOLVER_CALLBACK_HH__ */
