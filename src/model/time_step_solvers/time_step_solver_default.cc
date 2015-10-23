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
#include "sparse_matrix_aij.hh"

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
                                             const TimeStepSolverType & type,
                                             const ID & id, UInt memory_id)
    : TimeStepSolver(dof_manager, type, id, memory_id),
      dof_manager(dof_manager), is_mass_lumped(false) {

  switch (type) {
  case _tsst_forward_euler_lumped:
    this->is_mass_lumped = true;
  case _tsst_forward_euler: {

    this->integration_scheme = new ForwardEuler(dof_manager);
    break;
  }
  case _tsst_trapezoidal_rule_1: {
    this->integration_scheme = new TrapezoidalRule1(dof_manager);
    break;
  }
  case _tsst_backward_euler: {
    this->integration_scheme = new BackwardEuler(dof_manager);
    break;
  }
  case _tsst_central_difference_lumped:
    this->is_mass_lumped = true;
  case _tsst_central_difference: {
    this->integration_scheme = new CentralDifference(dof_manager);
    break;
  }
  case _tsst_fox_goodwin: {
    this->integration_scheme = new FoxGoodwin(dof_manager);
    break;
  }
  case _tsst_trapezoidal_rule_2: {
    this->integration_scheme = new TrapezoidalRule2(dof_manager);
    break;
  }
  case _tsst_linear_acceleration: {
    this->integration_scheme = new LinearAceleration(dof_manager);
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

  TimeStepSolver::predictor();

  this->integration_scheme->predictor(this->dof_id, this->time_step);

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

  this->integration_scheme->corrector(
      IntegrationScheme::SolutionType(this->solution_type), this->dof_id,
      this->time_step);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void TimeStepSolverDefault::assembleJacobian() {
  AKANTU_DEBUG_IN();

  TimeStepSolver::assembleJacobian();

  this->integration_scheme->assembleJacobian(
      IntegrationScheme::SolutionType(this->solution_type), this->time_step);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void TimeStepSolverDefault::assembleResidual() {
  AKANTU_DEBUG_IN();

  TimeStepSolver::assembleResidual();

  //   if (!this->is_mass_lumped) {

  //     Array<Real> * Ma = new Array<Real>(*acceleration, true, "Ma");
  //     *Ma *= *mass_matrix;
  //     /// \todo check unit conversion for implicit dynamics
  //     //      *Ma /= f_m2a
  //     *residual -= *Ma;
  //     delete Ma;
  //   } else if (mass) {

  //     // else lumped mass
  //     UInt nb_nodes = acceleration->getSize();
  //     UInt nb_degree_of_freedom = acceleration->getNbComponent();

  //     Real * mass_val = mass->storage();
  //     Real * accel_val = acceleration->storage();
  //     Real * res_val = residual->storage();
  //     bool * blocked_dofs_val = blocked_dofs->storage();

  //     for (UInt n = 0; n < nb_nodes * nb_degree_of_freedom; ++n) {
  //       if (!(*blocked_dofs_val)) {
  //         *res_val -= *accel_val * *mass_val / f_m2a;
  //       }
  //       blocked_dofs_val++;
  //       res_val++;
  //       mass_val++;
  //       accel_val++;
  //     }
  //   } else {
  //     AKANTU_DEBUG_ERROR("No function called to assemble the mass matrix.");
  //   }

  //   // f -= Cv
  //   if (velocity_damping_matrix) {
  //     Array<Real> * Cv = new Array<Real>(*velocity);
  //     *Cv *= *velocity_damping_matrix;
  //     *residual -= *Cv;
  //     delete Cv;
  //   }
  // }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__
