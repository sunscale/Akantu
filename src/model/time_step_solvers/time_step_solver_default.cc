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
#include "mesh.hh"
#include "dof_manager_default.hh"
#include "sparse_matrix_aij.hh"

#include "pseudo_time.hh"
#include "integration_scheme_1st_order.hh"
#include "integration_scheme_2nd_order.hh"
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
    /// initialize a static time solver for allback dofs
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
    // Write a c++11 version of the constructor with initializer list that
    // contains the arguments for the integration scheme
    case _ist_generalized_trapezoidal:
    case _ist_newmark_beta:
      AKANTU_EXCEPTION(
          "This time step solvers cannot be created with this constructor");
    }
  }

  this->integration_schemes[dof_id] = integration_scheme;
  this->solution_types[dof_id] = solution_type;

  this->integration_schemes_owner.insert(dof_id);
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
    integration_scheme_it->second->predictor(this->time_step);
  }

  // UInt nb_degree_of_freedom = u.getSize() * u.getNbComponent();

  // Array<Real>::scalar_iterator incr_it =
  //     increment.begin_reinterpret(nb_degree_of_freedom);
  // Array<Real>::const_scalar_iterator u_it =
  //     u.begin_reinterpret(nb_degree_of_freedom);
  // Array<Real>::const_scalar_iterator u_end =
  //     u.end_reinterpret(nb_degree_of_freedom);

  // for (; u_it != u_end; ++u_it, ++incr_it) {
  //   *incr_it = *u_it - *incr_it;
  // }

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

    integration_scheme_it->second->corrector(
        solution_type, this->time_step);
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
