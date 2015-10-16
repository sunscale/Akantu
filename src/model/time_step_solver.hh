/**
 * @file   time_step_solver.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Mon Aug 24 12:42:04 2015
 *
 * @brief  This corresponding to the time step evolution solver
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
#include "aka_memory.hh"
#include "aka_array.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_TIME_STEP_SOLVER_HH__
#define __AKANTU_TIME_STEP_SOLVER_HH__

namespace akantu {
  class DOFManager;
}

__BEGIN_AKANTU__

class TimeStepSolver : public Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  TimeStepSolver(DOFManager & dof_manager, const ID & dof_id,
                 const TimeStepSolverType & type, const ID & id,
                 UInt memory_id);
  virtual ~TimeStepSolver();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void solveStep();

private:
  virtual void predictor() = 0;

  virtual void corrector() = 0;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// Underlying dof manager containing the dof to treat
  DOFManager & _dof_manager;

  /// The dof to make evolve
  ID dof_id;

  /// Type of solver
  TimeStepSolverType type;

  /// The increment of the dof from @f[ u_{n+1} = u_{n} + inc @f]
  Array<Real> increment;

  /// The time step for this solver
  Real time_step;
};

__END_AKANTU__

#endif /* __AKANTU_TIME_STEP_SOLVER_HH__ */
