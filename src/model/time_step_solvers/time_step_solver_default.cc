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
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
// void TimeStepSolverDefault::updateAcceleration() {
//   AKANTU_DEBUG_IN();

//   updateResidualInternal();

//   if (method == _explicit_lumped_mass) {
//     /* residual = residual_{n+1} - M * acceleration_n therefore
//        solution = increment acceleration not acceleration */
//     solveLumped(*increment_acceleration, *mass, *residual, *blocked_dofs,
//                 f_m2a);
//   } else if (method == _explicit_consistent_mass) {
//     solve<NewmarkBeta::_acceleration_corrector>(*increment_acceleration);
//   }

//   AKANTU_DEBUG_OUT();
// }
/* -------------------------------------------------------------------------- */
TimeStepSolverDefault::TimeStepSolverDefault(DOFManagerDefault & dof_manager,
                                             const ID & dof_id,
                                             const TimeStepSolverType & type,
                                             const ID & id, UInt memory_id)
    : TimeStepSolver(dof_manager, dof_id, type, id, memory_id),
      dof_manager(dof_manager) {
  switch (type) {
  case _tsst_forward_euler: {
    this->integration_scheme = new ForwardEuler();
    break;
  }
  case _tsst_trapezoidal_rule_1: {
    this->integration_scheme = new TrapezoidalRule1();
    break;
  }
  case _tsst_backward_euler: {
    this->integration_scheme = new BackwardEuler();
    break;
  }
  case _tsst_central_difference: {
    this->integration_scheme = new CentralDifference();
    break;
  }
  case _tsst_fox_goodwin: {
    this->integration_scheme = new FoxGoodwin();
    break;
  }
  case _tsst_trapezoidal_rule_2: {
    this->integration_scheme = new TrapezoidalRule2();
    break;
  }
  case _tsst_linear_acceleration: {
    this->integration_scheme = new LinearAceleration();
    break;
  }
  // Write a c++11 version of the constructor with initializer list that
  // contains the arguments for the integration scheme
  case _tsst_generalized_trapezoidal:
  case _tsst_newmark_beta:
    AKANTU_EXCEPTION(
        "This time step solvers cannot be created with this constructor");
  }
}

/* -------------------------------------------------------------------------- */
TimeStepSolverDefault::~TimeStepSolverDefault() {
  delete this->integration_scheme;
}

/* -------------------------------------------------------------------------- */
void TimeStepSolverDefault::predictor() {
  AKANTU_DEBUG_IN();

  Array<Real> & u = this->dof_manager.getDOFs(this->dof_id);
  const Array<bool> & blocked_dofs =
      this->dof_manager.getBlockedDOFs(this->dof_id);

  // increment.copy(u);

  if (this->integration_scheme->getOrder() == 1) {
    Array<Real> & u_dot = dof_manager.getDOFsDerivatives(this->dof_id, 1);
    IntegrationScheme1stOrder & int_scheme =
        *dynamic_cast<IntegrationScheme1stOrder *>(this->integration_scheme);
    int_scheme.integrationSchemePred(this->time_step, u, u_dot, blocked_dofs);
  } else if (this->integration_scheme->getOrder() == 2) {

    Array<Real> & u_dot = dof_manager.getDOFsDerivatives(this->dof_id, 1);
    Array<Real> & u_dot_dot = dof_manager.getDOFsDerivatives(this->dof_id, 2);

    IntegrationScheme2ndOrder & int_scheme =
        *dynamic_cast<IntegrationScheme2ndOrder *>(this->integration_scheme);
    int_scheme.integrationSchemePred(this->time_step, u, u_dot, u_dot_dot,
                                     blocked_dofs);
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

  Array<Real> & u = this->dof_manager.getDOFs(this->dof_id);
  const Array<Real> & solution = this->dof_manager.getSolution(this->dof_id);
  const Array<bool> & blocked_dofs =
      this->dof_manager.getBlockedDOFs(this->dof_id);

  // increment.copy(u);

  if (this->integration_scheme->getOrder() == 1) {
    Array<Real> & u_dot = dof_manager.getDOFsDerivatives(this->dof_id, 1);
    IntegrationScheme1stOrder & int_scheme =
        *dynamic_cast<IntegrationScheme1stOrder *>(this->integration_scheme);

    switch (this->corrector_type) {
    case IntegrationScheme1stOrder::_temperature_corrector:
      int_scheme.integrationSchemeCorrTemp(this->time_step, u, u_dot,
                                           blocked_dofs, solution);
      break;
    case IntegrationScheme1stOrder::_temperature_rate_corrector:
      int_scheme.integrationSchemeCorrTempRate(this->time_step, u, u_dot,
                                               blocked_dofs, solution);
      break;
    }
  } else if (this->integration_scheme->getOrder() == 2) {

    Array<Real> & u_dot = dof_manager.getDOFsDerivatives(this->dof_id, 1);
    Array<Real> & u_dot_dot = dof_manager.getDOFsDerivatives(this->dof_id, 2);

    IntegrationScheme2ndOrder & int_scheme =
        *dynamic_cast<IntegrationScheme2ndOrder *>(this->integration_scheme);

    switch (this->corrector_type) {
    case IntegrationScheme2ndOrder::_displacement_corrector:
      int_scheme.integrationSchemeCorrDispl(this->time_step, u, u_dot,
                                            u_dot_dot, blocked_dofs, solution);
      break;
    case IntegrationScheme2ndOrder::_velocity_corrector:
      int_scheme.integrationSchemeCorrVeloc(this->time_step, u, u_dot,
                                            u_dot_dot, blocked_dofs, solution);
      break;
    case IntegrationScheme2ndOrder::_acceleration_corrector:
      int_scheme.integrationSchemeCorrAccel(this->time_step, u, u_dot,
                                            u_dot_dot, blocked_dofs, solution);
      break;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void TimeStepSolverDefault::solveStep() {
  AKANTU_DEBUG_IN();

  // EventManager::sendEvent(
  //     SolidMechanicsModelEvent::BeforeSolveStepEvent(method));

  // this->predictor();
  // this->solver->solve();
  // this->corrector();

  // EventManager::sendEvent(
  //     SolidMechanicsModelEvent::AfterSolveStepEvent(method));

  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
