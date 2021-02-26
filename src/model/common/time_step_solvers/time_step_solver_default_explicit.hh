/**
 * @file   time_step_solver_default_explicit.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  Default solver for explicit resolution
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

#ifndef AKANTU_TIME_STEP_SOLVER_DEFAULT_EXPLICIT_HH_
#define AKANTU_TIME_STEP_SOLVER_DEFAULT_EXPLICIT_HH_

namespace akantu {

class TimeStepSolverDefaultExplicit : public TimeStepSolverDefault {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  TimeStepSolverDefaultExplicit();
  virtual ~TimeStepSolverDefaultExplicit();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void solveStep();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */

public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream,
                                 const TimeStepSolverDefaultExplicit & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu


//#include "time_step_solver_default_explicit_inline_impl.hh"

#endif /* AKANTU_TIME_STEP_SOLVER_DEFAULT_EXPLICIT_HH_ */
