/**
 * @file   time_step_solver_default.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Feb 21 2018
 *
 * @brief  Default implementation for the time stepper
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
#include "integration_scheme.hh"
#include "time_step_solver.hh"
/* -------------------------------------------------------------------------- */
#include <map>
#include <set>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_TIME_STEP_SOLVER_DEFAULT_HH__
#define __AKANTU_TIME_STEP_SOLVER_DEFAULT_HH__

namespace akantu {
class DOFManager;
}

namespace akantu {

class TimeStepSolverDefault : public TimeStepSolver {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  TimeStepSolverDefault(DOFManager & dof_manager,
                        const TimeStepSolverType & type,
                        NonLinearSolver & non_linear_solver,
                        SolverCallback & solver_callback, const ID & id,
                        UInt memory_id);

  ~TimeStepSolverDefault() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  /// registers an integration scheme for a given dof
  void
  setIntegrationSchemeInternal(const ID & dof_id,
                               const IntegrationSchemeType & type,
                               IntegrationScheme::SolutionType solution_type =
                                   IntegrationScheme::_not_defined) override;

public:
  bool hasIntegrationScheme(const ID & dof_id) const override;

  /// implementation of the TimeStepSolver::predictor()
  void predictor() override;
  /// implementation of the TimeStepSolver::corrector()
  void corrector() override;
  /// implementation of the TimeStepSolver::assembleMatrix()
  void assembleMatrix(const ID & matrix_id) override;

  //  void assembleLumpedMatrix(const ID & matrix_id) override;
  /// implementation of the TimeStepSolver::assembleResidual()
  void assembleResidual() override;
  void assembleResidual(const ID & residual_part) override;

  void beforeSolveStep() override;
  void afterSolveStep(bool converged = true) override;

  /// implementation of the generic TimeStepSolver::solveStep()
  void solveStep(SolverCallback & solver_callback) override;

private:

  template<class Func>
  void for_each_integrator(Func && function) {
    for (auto & pair : this->integration_schemes) {
      auto & dof_id = pair.first;
      auto & integration_scheme = pair.second;
      function(dof_id, *integration_scheme);
    }
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  using DOFsIntegrationSchemes =
      std::map<ID, std::unique_ptr<IntegrationScheme>>;
  using DOFsIntegrationSchemesSolutionTypes =
      std::map<ID, IntegrationScheme::SolutionType>;
  using DOFsIntegrationSchemesOwner = std::set<ID>;

  /// Underlying integration scheme per dof, \todo check what happens in dynamic
  /// in case of coupled equations
  DOFsIntegrationSchemes integration_schemes;

  /// defines if the solver is owner of the memory or not
  DOFsIntegrationSchemesOwner integration_schemes_owner;

  /// Type of corrector to use
  DOFsIntegrationSchemesSolutionTypes solution_types;

  /// define if the mass matrix is lumped or not
  bool is_mass_lumped{false};
};

} // namespace akantu

#endif /* __AKANTU_TIME_STEP_SOLVER_DEFAULT_HH__ */
