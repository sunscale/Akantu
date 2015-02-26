/**
 * @file   solid_mechanics_model_cohesive_inline_impl.cc
 *
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 *
 * @date   Thu Feb 26 14:23:45 2015
 *
 * @brief  Implementation of inline functions for the Cohesive element model
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

/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_SOLID_MECHANICS_MODEL_COHESIVE_INLINE_IMPL_CC__
#define __AKANTU_SOLID_MECHANICS_MODEL_COHESIVE_INLINE_IMPL_CC__

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<SolveConvergenceMethod cmethod, SolveConvergenceCriteria criteria>
bool SolidMechanicsModel::solveStepCohesive(Real tolerance, Real & error, UInt max_iteration,
                                            bool do_not_factorize) {
  EventManager::sendEvent(SolidMechanicsModelEvent::BeforeSolveStepEvent(method));
  this->implicitPred();

  bool something_converged = false;

  Array<Real> * displacement_tmp = NULL;
  Array<Real> * velocity_tmp = NULL;
  Array<Real> * acceleration_tmp = NULL;

  while(something_converged) {
    if(extrinsic) {
      // If in extrinsic saves the current displacements, velocities and accelerations
      Array<Real> * tmp_swap;

      if(!displacement_tmp) {
        displacement_tmp = new Array<Real>(*(this->displacement));
      } else {
        displacement_tmp->resize(this->displacement->getSize());
        displacement_tmp->copy(*(this->displacement));
      }
      tmp_swap = *displacement_tmp;
      displacement_tmp = this->displacement;
      this->displacement = tmp_swap;

      if(!velocity_tmp) {
        velocity_tmp = new Array<Real>(*(this->velocity));
      } else {
        velocity_tmp->resize(this->velocity->getSize());
        velocity_tmp->copy(*(this->velocity));
      }
      tmp_swap = *velocity_tmp;
      velocity_tmp = this->velocity;
      this->velocity = tmp_swap;

      if(!acceleration_tmp) {
        acceleration_tmp = new Array<Real>(*(this->acceleration));
      } else {
        acceleration_tmp->resize(this->acceleration->getSize());
        acceleration_tmp->copy(*(this->acceleration));
      }
      tmp_swap = *acceleration_tmp;
      acceleration_tmp = this->acceleration;
      this->acceleration = tmp_swap;
    }


    this->updateResidual();

    AKANTU_DEBUG_ASSERT(stiffness_matrix != NULL,
                        "You should first initialize the implicit solver and assemble the stiffness matrix");

    bool need_factorize = !do_not_factorize;

    if (method==_implicit_dynamic) {
      AKANTU_DEBUG_ASSERT(mass_matrix != NULL,
                          "You should first initialize the implicit solver and assemble the mass matrix");
    }

    switch (cmethod) {
    case _scm_newton_raphson_tangent:
    case _scm_newton_raphson_tangent_not_computed:
      break;
    case _scm_newton_raphson_tangent_modified:
      this->assembleStiffnessMatrix();
      break;
    default:
      AKANTU_DEBUG_ERROR("The resolution method " << cmethod << " has not been implemented!");
    }

    UInt iter = 0;
    bool converged = false;
    error = 0.;
    if(criteria == _scc_residual) {
      converged = this->testConvergence<criteria> (tolerance, error);
      if(converged) return converged;
    }

    do {
      if (cmethod == _scm_newton_raphson_tangent)
        this->assembleStiffnessMatrix();

      solve<NewmarkBeta::_displacement_corrector> (*increment, 1., need_factorize);

      this->implicitCorr();

      this->updateResidual();

      converged = this->testConvergence<criteria> (tolerance, error);

      iter++;
      AKANTU_DEBUG_INFO("[" << criteria << "] Convergence iteration "
                        << std::setw(std::log10(max_iteration)) << iter
                        << ": error " << error << (converged ? " < " : " > ") << tolerance);

      switch (cmethod) {
      case _scm_newton_raphson_tangent:
        need_factorize = true;
        break;
      case _scm_newton_raphson_tangent_not_computed:
      case _scm_newton_raphson_tangent_modified:
        need_factorize = false;
        break;
      default:
        AKANTU_DEBUG_ERROR("The resolution method " << cmethod << " has not been implemented!");
      }


    } while (!converged && iter < max_iteration);


    if (converged) {
      EventManager::sendEvent(SolidMechanicsModelEvent::AfterSolveStepEvent(method));
    } else if(iter == max_iteration) {
      AKANTU_DEBUG_WARNING("[" << criteria << "] Convergence not reached after "
                           << std::setw(std::log10(max_iteration)) << iter <<
                           " iteration" << (iter == 1 ? "" : "s") << "!" << std::endl);
    }

    if(extrinsic) {
      Array<Real> * tmp_swap;

      tmp_swap = *displacement_tmp;
      displacement_tmp = this->displacement;
      this->displacement = tmp_swap;

      tmp_swap = *velocity_tmp;
      velocity_tmp = this->velocity;
      this->velocity = tmp_swap;

      tmp_swap = *acceleration_tmp;
      acceleration_tmp = this->acceleration;
      this->acceleration = tmp_swap;

      UInt nb_cohesive_elements = this->mesh.getNbElements(this->spatial_dimension, _not_ghost, _ek_cohesive) +
        this->mesh.getNbElements(this->spatial_dimension, _ghost, _ek_cohesive);
      this->checkCohesiveStress();

      UInt new_nb_cohesive_elements = this->mesh.getNbElements(this->spatial_dimension, _not_ghost, _ek_cohesive) +
        this->mesh.getNbElements(this->spatial_dimension, _ghost, _ek_cohesive);

      if (new_nb_cohesive_elements != nb_cohesive_elements) something_converged = false;
    } else {
      something_converged = true;
    }
  }

  if(extrinsic) {
    this->displacement->copy(*displacement_tmp);
    this->velocity    ->copy(*velocity_tmp);
    this->acceleration->copy(*acceleration_tmp);

    delete displacement_tmp;
    delete velocity_tmp;
    delete acceleration_tmp;
  }

  return something_converged;
}

__END_AKANTU__

#if defined (AKANTU_PARALLEL_COHESIVE_ELEMENT)
#  include "solid_mechanics_model_cohesive_parallel_inline_impl.cc"
#endif

#endif /* __AKANTU_SOLID_MECHANICS_MODEL_COHESIVE_INLINE_IMPL_CC__ */
