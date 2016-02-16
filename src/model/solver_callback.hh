/**
 * @file   solver_callback.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Tue Sep 15 22:45:27 2015
 *
 * @brief  Class defining the interface for non_linear_solver callbacks
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
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */


#ifndef __AKANTU_SOLVER_CALLBACK_HH__
#define __AKANTU_SOLVER_CALLBACK_HH__

__BEGIN_AKANTU__

class SolverCallback {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  SolverCallback() {}
  virtual ~SolverCallback() {}

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// callback to assemble the Jacobian Matrix
  virtual void assembleJacobian() { AKANTU_DEBUG_TO_IMPLEMENT(); }

  /// callback to assemble the residual (rhs)
  virtual void assembleResidual() { AKANTU_DEBUG_TO_IMPLEMENT(); }

  /* ------------------------------------------------------------------------ */
  /* Dynamic simulations part                                                 */
  /* ------------------------------------------------------------------------ */
  /// callback for the predictor (in case of dynamic simulation)
  virtual void predictor() { AKANTU_DEBUG_TO_IMPLEMENT(); }

  /// callback for the corrector (in case of dynamic simulation)
  virtual void corrector() { AKANTU_DEBUG_TO_IMPLEMENT(); }
};

__END_AKANTU__


#endif /* __AKANTU_SOLVER_CALLBACK_HH__ */
