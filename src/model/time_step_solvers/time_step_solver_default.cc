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
void TimeStepSolverDefault::predictor() {
  AKANTU_DEBUG_IN();

  Array<Real> & u = dof_manager.getDOFs(this->dof_id);
  Array<Real> & blocked_dofs = dof_manager.getBlockedDOFs(this->dof_id);

  IncrementsMap::itreator inc_it = increments.find(dof_id);
  if(inc_it != increments.end()) {
    inc_it->second->copy(u);
  }

  if (this->integrator->getOrder() == 1) {
    Array<Real> & u_dot = dof_manager.getDOFsDerivatives(this->dof_id, 1);
    IntergrationScheme1stOrder & int =
        *dynamic_cast<IntergrationScheme1stOrder *>(integrator);
    int.integrationSchemePred(this->time_step, u, u_dot, blocked_dofs);
  } else if (this->integrator->getOrder() == 2) {

    Array<Real> & u_dot = dof_manager.getDOFsDerivatives(this->dof_id, 1);
    Array<Real> & u_dot_dot = dof_manager.getDOFsDerivatives(this->dof_id, 2);

    IntergrationScheme2ndOrder & int =
        *dynamic_cast<IntergrationScheme2ndOrder *>(integrator);
    int.integrationSchemePred(this->time_step, u, u_dot, u_dot_dot,
                              blocked_dofs);
  }

  if (inc_it != increments.end()) {
    UInt nb_degree_of_freedom =
        u.getSize() * u.getNbComponent();

    Array<Real>::scalar_iterator incr_it = inc_it->second->begin_reinterpret(nb_degree_of_freedom);
    Array<Real>::const_scalar_iterator u_it = u.begin_reinterpret(nb_degree_of_freedom);
    Array<Real>::const_scalar_iterator u_end = u.end_reinterpret(nb_degree_of_freedom);

    for (; u_it != u_end; ++u_it, ++incr_it) {
      *incr_it = *u_it - *incr_it;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void TimeStepSolverDefault::corrector() {
  AKANTU_DEBUG_IN();

  Array<Real> & u = dof_manager.getDOFs(this->dof_id);
  Array<Real> & blocked_dofs = dof_manager.getBlockedDOFs(this->dof_id);

  integrator->integrationSchemeCorrAccel(time_step, *displacement, *velocity,
                                         *acceleration, *blocked_dofs,
                                         *increment_acceleration);

  if (previous_dofs)
    previous_dofs.copy(u);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void TimeStepSolverDefault::solveStep() {
  AKANTU_DEBUG_IN();

  EventManager::sendEvent(
      SolidMechanicsModelEvent::BeforeSolveStepEvent(method));

  this->predictor();
  this->solver->solve();
  this->corrector();

  EventManager::sendEvent(
      SolidMechanicsModelEvent::AfterSolveStepEvent(method));

  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
