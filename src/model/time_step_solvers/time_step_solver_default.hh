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
#include "integration_scheme.hh"
#include "time_step_solver.hh"
/* -------------------------------------------------------------------------- */
#include <map>
#include <set>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_TIME_STEP_SOLVER_DEFAULT_HH__
#define __AKANTU_TIME_STEP_SOLVER_DEFAULT_HH__

namespace akantu {
class DOFManagerDefault;
}

namespace akantu {

class TimeStepSolverDefault : public TimeStepSolver {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  TimeStepSolverDefault(DOFManagerDefault & dof_manager,
                        const TimeStepSolverType & type,
                        NonLinearSolver & non_linear_solver, const ID & id,
                        UInt memory_id);

  ~TimeStepSolverDefault() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// registers an integration scheme for a given dof
  void setIntegrationScheme(const ID & dof_id,
                            const IntegrationSchemeType & type,
                            IntegrationScheme::SolutionType solution_type =
                                IntegrationScheme::_not_defined) override;
  bool hasIntegrationScheme(const ID & dof_id) const override;

  /// implementation of the TimeStepSolver::predictor()
  void predictor() override;
  /// implementation of the TimeStepSolver::corrector()
  void corrector() override;
  /// implementation of the TimeStepSolver::assembleMatrix()
  void assembleMatrix(const ID & matrix_id) override;
  /// implementation of the TimeStepSolver::assembleResidual()
  void assembleResidual() override;

  /// implementation of the generic TimeStepSolver::solveStep()
  void solveStep(SolverCallback & solver_callback) override;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  using DOFsIntegrationSchemes = std::map<ID, std::unique_ptr<IntegrationScheme>>;
  using DOFsIntegrationSchemesSolutionTypes = std::map<ID, IntegrationScheme::SolutionType>;
  using DOFsIntegrationSchemesOwner = std::set<ID>;

  /// DOFManager with its real type
  DOFManagerDefault & dof_manager;

  /// Underlying integration scheme per dof, \todo check what happens in dynamic
  /// in case of coupled equations
  DOFsIntegrationSchemes integration_schemes;

  /// defines if the solver is owner of the memory or not
  DOFsIntegrationSchemesOwner integration_schemes_owner;

  /// Type of corrector to use
  DOFsIntegrationSchemesSolutionTypes solution_types;

  /// define if the mass matrix is lumped or not
  bool is_mass_lumped;
};

} // namespace akantu

#endif /* __AKANTU_TIME_STEP_SOLVER_DEFAULT_HH__ */
