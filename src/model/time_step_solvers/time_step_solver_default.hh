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
  class DOFManagerDefault;
}

__BEGIN_AKANTU__

class TimeStepSolverDefault : public TimeStepSolver {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  TimeStepSolverDefault(DOFManagerDefault & dof_manager,
                        const TimeStepSolverType & type, const ID & id,
                        UInt memory_id);

  virtual ~TimeStepSolverDefault();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// implementation of the TimeStepSolver::predictor()
  virtual void predictor();
  /// implementation of the TimeStepSolver::corrector()
  virtual void corrector();
  /// implementation of the TimeStepSolver::assembleJacobian()
  virtual void assembleJacobian();
  /// implementation of the TimeStepSolver::assembleResidual()
  virtual void assembleResidual();

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
  /// DOFManager with its real type
  DOFManagerDefault & dof_manager;

  /// Underlying integration scheme
  IntegrationScheme * integration_scheme;

  /// Type of corrector to use
  UInt solution_type;

  /// define if the mass matrix is lumped or not
  bool is_mass_lumped;
};

__END_AKANTU__

#endif /* __AKANTU_TIME_STEP_SOLVER_DEFAULT_HH__ */
