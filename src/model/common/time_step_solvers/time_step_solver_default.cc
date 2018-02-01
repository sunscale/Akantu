/**
 * @file   time_step_solver_default.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Sep 15 2015
 * @date last modification: Wed Feb 21 2018
 *
 * @brief  Default implementation of the time step solver
 *
 * @section LICENSE
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "time_step_solver_default.hh"
#include "dof_manager_default.hh"
#include "integration_scheme_1st_order.hh"
#include "integration_scheme_2nd_order.hh"
#include "mesh.hh"
#include "non_linear_solver.hh"
#include "pseudo_time.hh"
#include "sparse_matrix_aij.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
TimeStepSolverDefault::TimeStepSolverDefault(
    DOFManagerDefault & dof_manager, const TimeStepSolverType & type,
    NonLinearSolver & non_linear_solver, const ID & id, UInt memory_id)
    : TimeStepSolver(dof_manager, type, non_linear_solver, id, memory_id),
      dof_manager(dof_manager), is_mass_lumped(false) {
  switch (type) {
  case TimeStepSolverType::_dynamic:
    break;
  case TimeStepSolverType::_dynamic_lumped:
    this->is_mass_lumped = true;
    break;
  case TimeStepSolverType::_static:
    /// initialize a static time solver for callback dofs
    break;
  default:
    AKANTU_TO_IMPLEMENT();
  }
}

/* -------------------------------------------------------------------------- */
void TimeStepSolverDefault::setIntegrationScheme(
    const ID & dof_id, const IntegrationSchemeType & type,
    IntegrationScheme::SolutionType solution_type) {
  if (this->integration_schemes.find(dof_id) !=
      this->integration_schemes.end()) {
    AKANTU_EXCEPTION("Their DOFs "
                     << dof_id
                     << "  have already an integration scheme associated");
  }

  std::unique_ptr<IntegrationScheme> integration_scheme;
  if (this->is_mass_lumped) {
    switch (type) {
    case IntegrationSchemeType::_forward_euler: {
      integration_scheme = std::make_unique<ForwardEuler>(dof_manager, dof_id);
      break;
    }
    case IntegrationSchemeType::_central_difference: {
      integration_scheme =
          std::make_unique<CentralDifference>(dof_manager, dof_id);
      break;
    }
    default:
      AKANTU_EXCEPTION(
          "This integration scheme cannot be used in lumped dynamic");
    }
  } else {
    switch (type) {
    case IntegrationSchemeType::_pseudo_time: {
      integration_scheme = std::make_unique<PseudoTime>(dof_manager, dof_id);
      break;
    }
    case IntegrationSchemeType::_forward_euler: {
      integration_scheme = std::make_unique<ForwardEuler>(dof_manager, dof_id);
      break;
    }
    case IntegrationSchemeType::_trapezoidal_rule_1: {
      integration_scheme =
          std::make_unique<TrapezoidalRule1>(dof_manager, dof_id);
      break;
    }
    case IntegrationSchemeType::_backward_euler: {
      integration_scheme = std::make_unique<BackwardEuler>(dof_manager, dof_id);
      break;
    }
    case IntegrationSchemeType::_central_difference: {
      integration_scheme =
          std::make_unique<CentralDifference>(dof_manager, dof_id);
      break;
    }
    case IntegrationSchemeType::_fox_goodwin: {
      integration_scheme = std::make_unique<FoxGoodwin>(dof_manager, dof_id);
      break;
    }
    case IntegrationSchemeType::_trapezoidal_rule_2: {
      integration_scheme =
          std::make_unique<TrapezoidalRule2>(dof_manager, dof_id);
      break;
    }
    case IntegrationSchemeType::_linear_acceleration: {
      integration_scheme =
          std::make_unique<LinearAceleration>(dof_manager, dof_id);
      break;
    }
    case IntegrationSchemeType::_generalized_trapezoidal: {
      integration_scheme =
          std::make_unique<GeneralizedTrapezoidal>(dof_manager, dof_id);
      break;
    }
    case IntegrationSchemeType::_newmark_beta:
      integration_scheme = std::make_unique<NewmarkBeta>(dof_manager, dof_id);
      break;
    }
  }

  AKANTU_DEBUG_ASSERT(integration_scheme,
                      "No integration scheme was found for the provided types");

  auto && matrices_names = integration_scheme->getNeededMatrixList();
  for (auto && name : matrices_names) {
    needed_matrices.insert({name, _mt_not_defined});
  }

  this->integration_schemes[dof_id] = std::move(integration_scheme);
  this->solution_types[dof_id] = solution_type;

  this->integration_schemes_owner.insert(dof_id);
}

/* -------------------------------------------------------------------------- */
bool TimeStepSolverDefault::hasIntegrationScheme(const ID & dof_id) const {
  return this->integration_schemes.find(dof_id) !=
         this->integration_schemes.end();
}

/* -------------------------------------------------------------------------- */
TimeStepSolverDefault::~TimeStepSolverDefault() = default;

/* -------------------------------------------------------------------------- */
void TimeStepSolverDefault::solveStep(SolverCallback & solver_callback) {
  this->solver_callback = &solver_callback;

  this->solver_callback->beforeSolveStep();
  this->non_linear_solver.solve(*this);
  this->solver_callback->afterSolveStep();

  this->solver_callback = nullptr;
}

/* -------------------------------------------------------------------------- */
void TimeStepSolverDefault::predictor() {
  TimeStepSolver::predictor();

  for (auto && pair : this->integration_schemes) {
    const auto & dof_id = pair.first;
    auto & integration_scheme = pair.second;

    if (this->dof_manager.hasPreviousDOFs(dof_id)) {
      this->dof_manager.savePreviousDOFs(dof_id);
    }

    /// integrator predictor
    integration_scheme->predictor(this->time_step);
  }
}

/* -------------------------------------------------------------------------- */
void TimeStepSolverDefault::corrector() {
  AKANTU_DEBUG_IN();

  TimeStepSolver::corrector();

  for (auto & pair : this->integration_schemes) {
    auto & dof_id = pair.first;
    auto & integration_scheme = pair.second;

    const auto & solution_type = this->solution_types[dof_id];
    integration_scheme->corrector(solution_type, this->time_step);

    /// computing the increment of dof if needed
    if (this->dof_manager.hasDOFsIncrement(dof_id)) {
      if (!this->dof_manager.hasPreviousDOFs(dof_id)) {
        AKANTU_DEBUG_WARNING("In order to compute the increment of "
                             << dof_id << " a 'previous' has to be registered");
        continue;
      }

      Array<Real> & increment = this->dof_manager.getDOFsIncrement(dof_id);
      Array<Real> & previous = this->dof_manager.getPreviousDOFs(dof_id);

      UInt dof_array_comp = this->dof_manager.getDOFs(dof_id).getNbComponent();

      auto prev_dof_it = previous.begin(dof_array_comp);
      auto incr_it = increment.begin(dof_array_comp);
      auto incr_end = increment.end(dof_array_comp);

      increment.copy(this->dof_manager.getDOFs(dof_id));
      for (; incr_it != incr_end; ++incr_it, ++prev_dof_it) {
        *incr_it -= *prev_dof_it;
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void TimeStepSolverDefault::assembleMatrix(const ID & matrix_id) {
  AKANTU_DEBUG_IN();

  TimeStepSolver::assembleMatrix(matrix_id);

  if (matrix_id != "J")
    return;

  for (auto & pair : this->integration_schemes) {
    auto & dof_id = pair.first;
    auto & integration_scheme = pair.second;

    const auto & solution_type = this->solution_types[dof_id];

    integration_scheme->assembleJacobian(solution_type, this->time_step);
  }

  this->dof_manager.applyBoundary("J");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void TimeStepSolverDefault::assembleResidual() {
  if (this->needed_matrices.find("M") != needed_matrices.end()) {
    if (this->is_mass_lumped) {
      this->assembleLumpedMatrix("M");
    } else {
      this->assembleMatrix("M");
    }
  }

  TimeStepSolver::assembleResidual();

  for (auto && pair : this->integration_schemes) {
    auto & integration_scheme = pair.second;

    integration_scheme->assembleResidual(this->is_mass_lumped);
  }
}

/* -------------------------------------------------------------------------- */
void TimeStepSolverDefault::assembleResidual(const ID & residual_part) {
  AKANTU_DEBUG_IN();

  if (this->needed_matrices.find("M") != needed_matrices.end()) {
    if (this->is_mass_lumped) {
      this->assembleLumpedMatrix("M");
    } else {
      this->assembleMatrix("M");
    }
  }

  if (residual_part != "inertial") {
    TimeStepSolver::assembleResidual(residual_part);
  }

  if (residual_part == "inertial") {
    for (auto & pair : this->integration_schemes) {
      auto & integration_scheme = pair.second;

      integration_scheme->assembleResidual(this->is_mass_lumped);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
