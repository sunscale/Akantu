/**
 * @file   solid_mechanics_model_cohesive_inline_impl.cc
 *
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Jan 18 2013
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Implementation of inline functions for the Cohesive element model
 *
 * @section LICENSE
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_cohesive.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SOLID_MECHANICS_MODEL_COHESIVE_INLINE_IMPL_CC__
#define __AKANTU_SOLID_MECHANICS_MODEL_COHESIVE_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
// template <SolveConvergenceMethod cmethod, SolveConvergenceCriteria criteria>
// bool SolidMechanicsModelCohesive::solveStepCohesive(
//     Real tolerance, Real & error, UInt max_iteration, bool load_reduction,
//     Real tol_increase_factor, bool do_not_factorize) {

// // EventManager::sendEvent(
// //     SolidMechanicsModelEvent::BeforeSolveStepEvent(method));
// // this->implicitPred();

// // bool insertion_new_element = true;
// // bool converged = false;
// // Array<Real> * displacement_tmp = NULL;
// // Array<Real> * velocity_tmp = NULL;
// // Array<Real> * acceleration_tmp = NULL;
// // StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
// // Int prank = comm.whoAmI();

// // /// Loop for the insertion of new cohesive elements
// // while (insertion_new_element) {

// //   if (is_extrinsic) {
// //     /**
// //      * If in extrinsic the solution of the previous incremental step
// //      * is saved in temporary arrays created for displacements,
// //      * velocities and accelerations. Such arrays are used to find
// //      * the solution with the Newton-Raphson scheme (this is done by
// //      * pointing the pointer "displacement" to displacement_tmp). In
// //      * this way, inside the array "displacement" is kept the
// //      * solution of the previous incremental step, and in
// //      * "displacement_tmp" is saved the current solution.
// //      */

// //     if (!displacement_tmp)
// //       displacement_tmp = new Array<Real>(0, spatial_dimension);

// //     displacement_tmp->copy(*(this->displacement));

// //     if (!velocity_tmp)
// //       velocity_tmp = new Array<Real>(0, spatial_dimension);

// //     velocity_tmp->copy(*(this->velocity));

// //     if (!acceleration_tmp) {
// //       acceleration_tmp = new Array<Real>(0, spatial_dimension);
// //     }

// //     acceleration_tmp->copy(*(this->acceleration));

// //     std::swap(displacement, displacement_tmp);
// //     std::swap(velocity, velocity_tmp);
// //     std::swap(acceleration, acceleration_tmp);
// //   }

// //   this->updateResidual();

// //   AKANTU_DEBUG_ASSERT(stiffness_matrix != NULL,
// //                       "You should first initialize the implicit solver and
// "
// //                       "assemble the stiffness matrix");

// //   bool need_factorize = !do_not_factorize;

// //   if (method == _implicit_dynamic) {
// //     AKANTU_DEBUG_ASSERT(mass_matrix != NULL, "You should first initialize
// "
// //                                              "the implicit solver and "
// //                                              "assemble the mass matrix");
// //   }

// //   switch (cmethod) {
// //   case _scm_newton_raphson_tangent:
// //   case _scm_newton_raphson_tangent_not_computed:
// //     break;
// //   case _scm_newton_raphson_tangent_modified:
// //     this->assembleStiffnessMatrix();
// //     break;
// //   default:
// //     AKANTU_ERROR("The resolution method "
// //                        << cmethod << " has not been implemented!");
// //   }

// //   UInt iter = 0;
// //   converged = false;
// //   error = 0.;
// //   if (criteria == SolveConvergenceCriteria::_residual) {
// //     converged = this->testConvergence<criteria>(tolerance, error);
// //     if (converged)
// //       return converged;
// //   }

// //   /// Loop to solve the nonlinear system
// //   do {
// //     if (cmethod == _scm_newton_raphson_tangent)
// //       this->assembleStiffnessMatrix();

// //     solve<NewmarkBeta::_displacement_corrector>(*increment, 1.,
// //                                                 need_factorize);

// //     this->implicitCorr();

// //     this->updateResidual();

// //     converged = this->testConvergence<criteria>(tolerance, error);

// //     iter++;
// //     AKANTU_DEBUG_INFO("[" << criteria << "] Convergence iteration "
// //                           << std::setw(std::log10(max_iteration)) << iter
// //                           << ": error " << error
// //                           << (converged ? " < " : " > ") << tolerance);

// //     switch (cmethod) {
// //     case _scm_newton_raphson_tangent:
// //       need_factorize = true;
// //       break;
// //     case _scm_newton_raphson_tangent_not_computed:
// //     case _scm_newton_raphson_tangent_modified:
// //       need_factorize = false;
// //       break;
// //     default:
// //       AKANTU_ERROR("The resolution method "
// //                          << cmethod << " has not been implemented!");
// //     }

// //   } while (!converged && iter < max_iteration);

// //   /**
// //    * This is to save the obtained result and proceed with the
// //    * simulation even if the error is higher than the pre-fixed
// //    * tolerance. This is done only after loading reduction
// //    * (load_reduction = true).
// //    */
// //   //    if (load_reduction && (error < tolerance * tol_increase_factor))
// //   //    converged = true;
// //   if ((error < tolerance * tol_increase_factor))
// //     converged = true;

// //   if (converged) {

// //   } else if (iter == max_iteration) {
// //     if (prank == 0) {
// //       AKANTU_DEBUG_WARNING(
// //           "[" << criteria << "] Convergence not reached after "
// //               << std::setw(std::log10(max_iteration)) << iter << "
// iteration"
// //               << (iter == 1 ? "" : "s") << "!" << std::endl);
// //     }
// //   }

// //   if (is_extrinsic) {
// //     /**
// //      * If is extrinsic the pointer "displacement" is moved back to
// //      * the array displacement. In this way, the array displacement is
// //      * correctly resized during the checkCohesiveStress function (in
// //      * case new cohesive elements are added). This is possible
// //      * because the procedure called by checkCohesiveStress does not
// //      * use the displacement field (the correct one is now stored in
// //      * displacement_tmp), but directly the stress field that is
// //      * already computed.
// //      */
// //     Array<Real> * tmp_swap;

// //     tmp_swap = displacement_tmp;
// //     displacement_tmp = this->displacement;
// //     this->displacement = tmp_swap;

// //     tmp_swap = velocity_tmp;
// //     velocity_tmp = this->velocity;
// //     this->velocity = tmp_swap;

// //     tmp_swap = acceleration_tmp;
// //     acceleration_tmp = this->acceleration;
// //     this->acceleration = tmp_swap;

// //     /// If convergence is reached, call checkCohesiveStress in order
// //     /// to check if cohesive elements have to be introduced
// //     if (converged) {

// //       UInt new_cohesive_elements = checkCohesiveStress();

// //       if (new_cohesive_elements == 0) {
// //         insertion_new_element = false;
// //       } else {
// //         insertion_new_element = true;
// //       }
// //     }
// //   }

// //   if (!converged && load_reduction)
// //     insertion_new_element = false;

// //   /**
// //    * If convergence is not reached, there is the possibility to
// //    * return back to the main file and reduce the load. Before doing
// //    * this, a pre-fixed value as to be defined for the parameter
// //    * delta_max of the cohesive elements introduced in the current
// //    * incremental step. This is done by calling the function
// //    * checkDeltaMax.
// //    */
// //   if (!converged) {
// //     insertion_new_element = false;

// //     for (UInt m = 0; m < materials.size(); ++m) {
// //       try {
// //         MaterialCohesive & mat =
// //             aka::as_type<MaterialCohesive>(*materials[m]);
// //         mat.checkDeltaMax(_not_ghost);
// //       } catch (std::bad_cast &) {
// //       }
// //     }
// //   }

// // } // end loop for the insertion of new cohesive elements

// // /**
// //  * When the solution to the current incremental step is computed (no
// //  * more cohesive elements have to be introduced), call the function
// //  * to compute the energies.
// //  */
// // if ((is_extrinsic && converged)) {

// //   for (UInt m = 0; m < materials.size(); ++m) {
// //     try {
// //       MaterialCohesive & mat =
// //           aka::as_type<MaterialCohesive>(*materials[m]);
// //       mat.computeEnergies();
// //     } catch (std::bad_cast & bce) {
// //     }
// //   }

// //   EventManager::sendEvent(
// //       SolidMechanicsModelEvent::AfterSolveStepEvent(method));

// //   /**
// //    * The function resetVariables is necessary to correctly set a
// //    * variable that permit to decrease locally the penalty parameter
// //    * for compression.
// //    */
// //   for (UInt m = 0; m < materials.size(); ++m) {
// //     try {
// //       MaterialCohesive & mat =
// //           aka::as_type<MaterialCohesive>(*materials[m]);
// //       mat.resetVariables(_not_ghost);
// //     } catch (std::bad_cast &) {
// //     }
// //   }

// //   /// The correct solution is saved
// //   this->displacement->copy(*displacement_tmp);
// //   this->velocity->copy(*velocity_tmp);
// //   this->acceleration->copy(*acceleration_tmp);
// // }

// // delete displacement_tmp;
// // delete velocity_tmp;
// // delete acceleration_tmp;

// // return insertion_new_element;
//}

} // namespace akantu

#endif /* __AKANTU_SOLID_MECHANICS_MODEL_COHESIVE_INLINE_IMPL_CC__ */
