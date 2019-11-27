/**
 * @file   time_step_solver.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Feb 21 2018
 *
 * @brief  This corresponding to the time step evolution solver
 *
 * @section LICENSE
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
#include "aka_array.hh"
#include "aka_memory.hh"
#include "integration_scheme.hh"
#include "parameter_registry.hh"
#include "solver_callback.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_TIME_STEP_SOLVER_HH__
#define __AKANTU_TIME_STEP_SOLVER_HH__

namespace akantu {
class DOFManager;
class NonLinearSolver;
} // namespace akantu

namespace akantu {

class TimeStepSolver : public Memory,
                       public ParameterRegistry,
                       public SolverCallback {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  TimeStepSolver(DOFManager & dof_manager, const TimeStepSolverType & type,
                 NonLinearSolver & non_linear_solver,
                 SolverCallback & solver_callback, const ID & id,
                 UInt memory_id);
  ~TimeStepSolver() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// solves on step
  virtual void solveStep(SolverCallback & solver_callback) = 0;

  /// register an integration scheme for a given dof
  void setIntegrationScheme(const ID & dof_id,
                            const IntegrationSchemeType & type,
                            IntegrationScheme::SolutionType solution_type =
                                IntegrationScheme::_not_defined);

protected:
  /// register an integration scheme for a given dof
  virtual void
  setIntegrationSchemeInternal(const ID & dof_id,
                               const IntegrationSchemeType & type,
                               IntegrationScheme::SolutionType solution_type =
                                   IntegrationScheme::_not_defined) = 0;

public:
  /// replies if a integration scheme has been set
  virtual bool hasIntegrationScheme(const ID & dof_id) const = 0;
  /* ------------------------------------------------------------------------ */
  /* Solver Callback interface                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// implementation of the SolverCallback::getMatrixType()
  MatrixType getMatrixType(const ID &) final { return _mt_not_defined; }
  /// implementation of the SolverCallback::predictor()
  void predictor() override;
  /// implementation of the SolverCallback::corrector()
  void corrector() override;
  /// implementation of the SolverCallback::assembleJacobian()
  void assembleMatrix(const ID & matrix_id) override;
  /// implementation of the SolverCallback::assembleJacobian()
  void assembleLumpedMatrix(const ID & matrix_id) override;
  /// implementation of the SolverCallback::assembleResidual()
  void assembleResidual() override;
  /// implementation of the SolverCallback::assembleResidual()
  void assembleResidual(const ID & residual_part) override;

  void beforeSolveStep() override;
  void afterSolveStep(bool converged = true) override;

  bool canSplitResidual() override {
    return solver_callback->canSplitResidual();
  }
  /* ------------------------------------------------------------------------ */
  /* Accessor                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(TimeStep, time_step, Real);
  AKANTU_SET_MACRO(TimeStep, time_step, Real);

  AKANTU_GET_MACRO(NonLinearSolver, non_linear_solver, const NonLinearSolver &);
  AKANTU_GET_MACRO_NOT_CONST(NonLinearSolver, non_linear_solver,
                             NonLinearSolver &);

protected:
  MatrixType getCommonMatrixType();

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// Underlying dof manager containing the dof to treat
  DOFManager & _dof_manager;

  /// Type of solver
  TimeStepSolverType type;

  /// The time step for this solver
  Real time_step;

  /// Temporary storage for solver callback
  SolverCallback * solver_callback;

  /// NonLinearSolver used by this tome step solver
  NonLinearSolver & non_linear_solver;

  /// List of required matrices
  std::map<std::string, MatrixType> needed_matrices;

  /// specifies if the solvers gives to full solution or just the increment of
  /// solution
  bool is_solution_increment{true};
};

} // namespace akantu

#endif /* __AKANTU_TIME_STEP_SOLVER_HH__ */
