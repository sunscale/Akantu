/**
 * @file   time_step_solver_default.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Mon Aug 24 17:10:29 2015
 *
 * @brief  Default implementation for the time stepper
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
#include "time_step_solver.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_TIME_STEP_SOLVER_DEFAULT_HH__
#define __AKANTU_TIME_STEP_SOLVER_DEFAULT_HH__

namespace akantu {
  class IntegrationScheme;
}

__BEGIN_AKANTU__

class TimeStepSolverDefault : public TimeStepSolver {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  TimeStepSolverDefault(const TimeStepSolverType & type);
  virtual ~TimeStepSolverDefault();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// implementation of the TimeStepSolver::predictor()
  virtual void predictor();
  /// implementation of the TimeStepSolver::corrector()
  virtual void corrector();

  /// implementation of the generic TimeStepSolver::solveStep()
  virtual void solveStep();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// Underlying integration scheme
  IntegrationScheme * integration_scheme;

  /// Type of corrector to use
  UInt corrector_type;
};

__END_AKANTU__


#endif /* __AKANTU_TIME_STEP_SOLVER_DEFAULT_HH__ */
