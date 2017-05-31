/**
 * @file   time_step_solver_default.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Sep 16 10:20:55 2015
 *
 * @brief  Default implementation of the time step solver
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
#include "time_step_solver_default.hh"
#include "dof_manager_default.hh"
#include "integration_scheme_1st_order.hh"
#include "integration_scheme_2nd_order.hh"
#include "mesh.hh"
#include "non_linear_solver.hh"
#include "pseudo_time.hh"
#include "sparse_matrix_aij.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
TimeStepSolverDefault::TimeStepSolverDefault(
    DOFManagerDefault & dof_manager, const TimeStepSolverType & type,
    NonLinearSolver & non_linear_solver, const ID & id, UInt memory_id)
    : TimeStepSolver(dof_manager, type, non_linear_solver, id, memory_id),
      dof_manager(dof_manager), is_mass_lumped(false) {
  switch (type) {
  case _tsst_dynamic:
    break;
  case _tsst_dynamic_lumped:
    this->is_mass_lumped = true;
    break;
  case _tsst_static:
    /// initialize a static time solver for callback dofs
    break;
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

  IntegrationScheme * integration_scheme = NULL;
  if (this->is_mass_lumped) {
    switch (type) {
    case _ist_forward_euler: {
      integration_scheme = new ForwardEuler(dof_manager, dof_id);
      break;
    }
    case _ist_central_difference: {
      integration_scheme = new CentralDifference(dof_manager, dof_id);
      break;
    }
    default:
      AKANTU_EXCEPTION(
          "This integration scheme cannot be used in lumped dynamic");
    }
  } else {
    switch (type) {
    case _ist_pseudo_time: {
      integration_scheme = new PseudoTime(dof_manager, dof_id);
      break;
    }
    case _ist_forward_euler: {
      integration_scheme = new ForwardEuler(dof_manager, dof_id);
      break;
    }
    case _ist_trapezoidal_rule_1: {
      integration_scheme = new TrapezoidalRule1(dof_manager, dof_id);
      break;
    }
    case _ist_backward_euler: {
      integration_scheme = new BackwardEuler(dof_manager, dof_id);
      break;
    }
    case _ist_central_difference: {
      integration_scheme = new CentralDifference(dof_manager, dof_id);
      break;
    }
    case _ist_fox_goodwin: {
      integration_scheme = new FoxGoodwin(dof_manager, dof_id);
      break;
    }
    case _ist_trapezoidal_rule_2: {
      integration_scheme = new TrapezoidalRule2(dof_manager, dof_id);
      break;
    }
    case _ist_linear_acceleration: {
      integration_scheme = new LinearAceleration(dof_manager, dof_id);
      break;
    }
    case _ist_generalized_trapezoidal: {
      integration_scheme = new GeneralizedTrapezoidal(dof_manager, dof_id);
      break;
    }
    case _ist_newmark_beta:
      integration_scheme = new NewmarkBeta(dof_manager, dof_id);
      break;
    }
  }

  AKANTU_DEBUG_ASSERT(integration_scheme != nullptr,
                      "No integration scheme was found for the provided types");
  this->integration_schemes[dof_id] = integration_scheme;
  this->solution_types[dof_id] = solution_type;

  this->integration_schemes_owner.insert(dof_id);
}

/* -------------------------------------------------------------------------- */
bool TimeStepSolverDefault::hasIntegrationScheme(const ID & dof_id) const {
  return this->integration_schemes.find(dof_id) !=
         this->integration_schemes.end();
}
/* -------------------------------------------------------------------------- */
TimeStepSolverDefault::~TimeStepSolverDefault() {
  DOFsIntegrationSchemesOwner::iterator it =
      this->integration_schemes_owner.begin();
  DOFsIntegrationSchemesOwner::iterator end =
      this->integration_schemes_owner.end();

  for (; it != end; ++it) {
    delete this->integration_schemes[*it];
  }
  this->integration_schemes.clear();
}

/* -------------------------------------------------------------------------- */
void TimeStepSolverDefault::solveStep(SolverCallback & solver_callback) {
  this->solver_callback = &solver_callback;

  this->non_linear_solver.solve(*this);

  this->solver_callback = NULL;
}

/* -------------------------------------------------------------------------- */
void TimeStepSolverDefault::predictor() {
  AKANTU_DEBUG_IN();

  TimeStepSolver::predictor();

  DOFsIntegrationSchemes::iterator integration_scheme_it =
      this->integration_schemes.begin();
  DOFsIntegrationSchemes::iterator integration_scheme_end =
      this->integration_schemes.end();

  for (; integration_scheme_it != integration_scheme_end;
       ++integration_scheme_it) {
    ID dof_id = integration_scheme_it->first;
    Array<Real> * previous = NULL;

    UInt dof_array_comp = this->dof_manager.getDOFs(dof_id).getNbComponent();

    if (this->dof_manager.hasPreviousDOFs(dof_id)) {
      this->dof_manager.savePreviousDOFs(dof_id);
    } else {
      if (this->dof_manager.hasDOFsIncrement(dof_id)) {
        previous = new Array<Real>(this->dof_manager.getDOFs(dof_id));
      }
    }

    /// integrator predictor
    integration_scheme_it->second->predictor(this->time_step);

    /// computing the increment of dof if needed
    if (this->dof_manager.hasDOFsIncrement(dof_id)) {
      Array<Real> & increment = this->dof_manager.getDOFsIncrement(dof_id);
      Array<Real>::vector_iterator incr_it = increment.begin(dof_array_comp);
      Array<Real>::vector_iterator incr_end = increment.end(dof_array_comp);
      Array<Real>::const_vector_iterator dof_it =
          this->dof_manager.getDOFs(dof_id).begin(dof_array_comp);
      Array<Real>::const_vector_iterator prev_dof_it;

      if (this->dof_manager.hasPreviousDOFs(dof_id)) {
        prev_dof_it =
            this->dof_manager.getPreviousDOFs(dof_id).begin(dof_array_comp);
      } else {
        prev_dof_it = previous->begin(dof_array_comp);
      }

      for (; incr_it != incr_end; ++incr_it) {
        *incr_it = *dof_it;
        *incr_it -= *prev_dof_it;
      }

      delete previous;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void TimeStepSolverDefault::corrector() {
  AKANTU_DEBUG_IN();

  TimeStepSolver::corrector();

  DOFsIntegrationSchemes::iterator integration_scheme_it =
      this->integration_schemes.begin();
  DOFsIntegrationSchemes::iterator integration_scheme_end =
      this->integration_schemes.end();

  for (; integration_scheme_it != integration_scheme_end;
       ++integration_scheme_it) {
    IntegrationScheme::SolutionType solution_type =
        this->solution_types[integration_scheme_it->first];

    integration_scheme_it->second->corrector(solution_type, this->time_step);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void TimeStepSolverDefault::assembleJacobian() {
  AKANTU_DEBUG_IN();

  TimeStepSolver::assembleJacobian();

  DOFsIntegrationSchemes::iterator integration_scheme_it =
      this->integration_schemes.begin();
  DOFsIntegrationSchemes::iterator integration_scheme_end =
      this->integration_schemes.end();

  for (; integration_scheme_it != integration_scheme_end;
       ++integration_scheme_it) {
    IntegrationScheme::SolutionType solution_type =
        this->solution_types[integration_scheme_it->first];

    integration_scheme_it->second->assembleJacobian(solution_type,
                                                    this->time_step);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void TimeStepSolverDefault::assembleResidual() {
  AKANTU_DEBUG_IN();

  TimeStepSolver::assembleResidual();

  DOFsIntegrationSchemes::iterator integration_scheme_it =
      this->integration_schemes.begin();
  DOFsIntegrationSchemes::iterator integration_scheme_end =
      this->integration_schemes.end();

  for (; integration_scheme_it != integration_scheme_end;
       ++integration_scheme_it) {
    integration_scheme_it->second->assembleResidual(this->is_mass_lumped);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__
