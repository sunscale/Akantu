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

  /// solve a step using a given pre instantiated time step solver and
  /// nondynamic linear solver
  void solveStep(const ID & time_step_solver_id,
                 const ID & non_linear_solver_id);

  /// Initialize a time solver that can be used afterwards with its id
  void initTimeStepSolver(const ID & time_step_solver_id, const ID & dof_id,
                          const TimeStepSolverType & time_step_solver_type);

  /* ------------------------------------------------------------------------ */
  /* NonLinearSolverCallback interface                                        */
  /* ------------------------------------------------------------------------ */
public:
  /// Predictor interface for the callback
  void predictor();

  /// Corrector interface for the callback
  void corrector();

  /// AssembleResidual interface for the callback
  virtual void assembleResidual() = 0;

  /// AssembleJacobian interface for the callback
  virtual void assembleJacobian() = 0;

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
  ID default_time_step_solver_id;

  /// List of instantiated non_linear_solvers
  std::set<ID> non_linear_solvers;

  /// Default nondynamic linear solver
  ID default_non_linear_solver_id;
};

__END_AKANTU__

#endif /* __AKANTU_MODEL_SOLVER_HH__ */
