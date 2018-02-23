/**
 * @file   non_linear_solver_callback.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Mon Sep 28 18:48:21 2015
 *
 * @brief  Interface to implement for the non linear solver to work
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

#ifndef __AKANTU_NON_LINEAR_SOLVER_CALLBACK_HH__
#define __AKANTU_NON_LINEAR_SOLVER_CALLBACK_HH__

namespace akantu {

class NonLinearSolverCallback {
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// callback to assemble the Jacobian Matrix
  virtual void assembleJacobian() { AKANTU_TO_IMPLEMENT(); }

  /// callback to assemble the residual (rhs)
  virtual void assembleResidual() { AKANTU_TO_IMPLEMENT(); }

  /* ------------------------------------------------------------------------ */
  /* Dynamic simulations part                                                 */
  /* ------------------------------------------------------------------------ */
  /// callback for the predictor (in case of dynamic simulation)
  virtual void predictor() { AKANTU_TO_IMPLEMENT(); }

  /// callback for the corrector (in case of dynamic simulation)
  virtual void corrector() { AKANTU_TO_IMPLEMENT(); }
};

} // akantu

#endif /* __AKANTU_NON_LINEAR_SOLVER_CALLBACK_HH__ */
