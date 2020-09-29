/**
 * @file   solver_callback.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Feb 21 2018
 *
 * @brief  Class defining the interface for non_linear_solver callbacks
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_SOLVER_CALLBACK_HH_
#define AKANTU_SOLVER_CALLBACK_HH_

namespace akantu {
class DOFManager;
}

namespace akantu {

class SolverCallback {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  explicit SolverCallback(DOFManager & dof_manager);
  explicit SolverCallback();
  /* ------------------------------------------------------------------------ */
  virtual ~SolverCallback();

protected:
  void setDOFManager(DOFManager & dof_manager);
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// get the type of matrix needed
  virtual MatrixType getMatrixType(const ID &) = 0;

  /// callback to assemble a Matrix
  virtual void assembleMatrix(const ID &) = 0;

  /// callback to assemble a lumped Matrix
  virtual void assembleLumpedMatrix(const ID &) = 0;

  /// callback to assemble the residual (rhs)
  virtual void assembleResidual() = 0;

  /// callback to assemble the rhs parts, (e.g. internal_forces +
  /// external_forces)
  virtual void assembleResidual(const ID & /*residual_part*/) {}

  /* ------------------------------------------------------------------------ */
  /* Dynamic simulations part                                                 */
  /* ------------------------------------------------------------------------ */
  /// callback for the predictor (in case of dynamic simulation)
  virtual void predictor() {}

  /// callback for the corrector (in case of dynamic simulation)
  virtual void corrector() {}

  /// tells if the residual can be computed in separated parts
  virtual bool canSplitResidual() { return false; }

  /* ------------------------------------------------------------------------ */
  /* management callbacks                                                     */
  /* ------------------------------------------------------------------------ */
  virtual void beforeSolveStep(){};
  virtual void afterSolveStep(bool /*converged*/ = true){};

protected:
  /// DOFManager prefixed to avoid collision in multiple inheritance cases
  DOFManager * sc_dof_manager{nullptr};
};

namespace debug {
  class SolverCallbackResidualPartUnknown : public Exception {
  public:
    SolverCallbackResidualPartUnknown(const ID & residual_part)
        : Exception(residual_part + " is not known here.") {}
  };
} // namespace debug

} // namespace akantu

#endif /* AKANTU_SOLVER_CALLBACK_HH_ */
