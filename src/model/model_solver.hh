/**
 * @file   model_solver.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Jul 22 10:53:10 2015
 *
 * @brief  Class regrouping the common solve interface to the different models
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
#include "parsable.hh"
#include "dof_manager.hh"
#include "non_linear_solver_callback.hh"
/* -------------------------------------------------------------------------- */
#include <set>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MODEL_SOLVER_HH__
#define __AKANTU_MODEL_SOLVER_HH__

__BEGIN_AKANTU__

class ModelSolver : public Parsable, public NonLinearSolverCallback {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  ModelSolver(const Mesh & mesh, const ID & id, UInt memory_id);
  virtual ~ModelSolver();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// Solve a step, implicit, explicit, static, ... any things based on the
  /// current configuration
  void solveStep();

  /// Initialize a time solver that can be used afterwards with its id
  void initTimeStepSolver(const ID & time_step_solver_id, const ID & dof_id,
                          const TimeStepSolverType & time_step_solver_type);

protected:
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
protected:
  DOFManager & getDOFManager() { return *this->dof_manager; }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// Underlying dof_manager (the brain...)
  DOFManager * dof_manager;

  /// List of instantiated time step solvers
  std::set<ID> time_step_solvers;

  /// Default time step solver to use
  ID default_time_step_solver;
};

__END_AKANTU__

#endif /* __AKANTU_MODEL_SOLVER_HH__ */
