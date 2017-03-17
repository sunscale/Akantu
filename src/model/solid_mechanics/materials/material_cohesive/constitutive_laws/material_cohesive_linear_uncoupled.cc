/**
 * @file   material_cohesive_linear_uncoupled.cc
 *
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 *
 * @date creation: Wed Feb 22 2012
 * @date last modification: Thu Jan 14 2016
 *
 * @brief  Linear irreversible cohesive law of mixed mode loading with
 * random stress definition for extrinsic type
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include <algorithm>
#include <numeric>

/* -------------------------------------------------------------------------- */
#include "material_cohesive_linear_uncoupled.hh"
#include "solid_mechanics_model_cohesive.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialCohesiveLinearUncoupled<spatial_dimension>
::MaterialCohesiveLinearUncoupled(SolidMechanicsModel & model,
				  const ID & id) :
  MaterialCohesiveLinear<spatial_dimension>(model,id),
  delta_n_max("delta_n_max", *this),
  delta_t_max("delta_t_max", *this),
  damage_n("damage_n", *this),
  damage_t("damage_t", *this){

  AKANTU_DEBUG_IN();

  this->registerParam("roughness", R, Real(1.),
                      _pat_parsable | _pat_readable,
                      "Roughness to define coupling between mode II and mode I");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveLinearUncoupled<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();

  MaterialCohesiveLinear<spatial_dimension>::initMaterial();

  delta_n_max.initialize(1);
  delta_t_max.initialize(1);
  damage_n.initialize(1);
  damage_t.initialize(1);
  
  delta_n_max.initializeHistory();
  delta_t_max.initializeHistory();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveLinearUncoupled<spatial_dimension>::computeTraction(const Array<Real> & normal,
                                                                ElementType el_type,
                                                                GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  delta_n_max.resize();
  delta_t_max.resize();
  damage_n.resize();
  damage_t.resize();
  
  /// define iterators
  Array<Real>::vector_iterator traction_it =
    this->tractions(el_type, ghost_type).begin(spatial_dimension);

  Array<Real>::vector_iterator traction_end =
    this->tractions(el_type, ghost_type).end(spatial_dimension);

  Array<Real>::vector_iterator opening_it =
    this->opening(el_type, ghost_type).begin(spatial_dimension);

  /// opening_prec is the opening of the previous step in the
  /// Newton-Raphson loop
  Array<Real>::vector_iterator opening_prec_it =
   this->opening_prec(el_type, ghost_type).begin(spatial_dimension);

  Array<Real>::vector_iterator contact_traction_it =
    this->contact_tractions(el_type, ghost_type).begin(spatial_dimension);

  Array<Real>::vector_iterator contact_opening_it =
    this->contact_opening(el_type, ghost_type).begin(spatial_dimension);

  Array<Real>::const_vector_iterator normal_it =
    this->normal.begin(spatial_dimension);

  Array<Real>::scalar_iterator sigma_c_it =
    this->sigma_c_eff(el_type, ghost_type).begin();

  Array<Real>::scalar_iterator delta_n_max_it =
    this->delta_n_max(el_type, ghost_type).begin();

  Array<Real>::scalar_iterator delta_t_max_it =
    this->delta_t_max(el_type, ghost_type).begin();

  Array<Real>::scalar_iterator delta_c_it =
    this->delta_c_eff(el_type, ghost_type).begin();

  Array<Real>::scalar_iterator damage_n_it =
    this->damage_n(el_type, ghost_type).begin();

  Array<Real>::scalar_iterator damage_t_it =
    this->damage_t(el_type, ghost_type).begin();

  Array<Real>::vector_iterator insertion_stress_it =
    this->insertion_stress(el_type, ghost_type).begin(spatial_dimension);

  Array<bool>::scalar_iterator reduction_penalty_it =
    this->reduction_penalty(el_type, ghost_type).begin();

  Vector<Real> normal_opening(spatial_dimension);
  Vector<Real> tangential_opening(spatial_dimension);

  if (!this->model->isExplicit()){
    this->delta_n_max(el_type, ghost_type).copy(this->delta_n_max.previous(el_type, ghost_type));
    this->delta_t_max(el_type, ghost_type).copy(this->delta_t_max.previous(el_type, ghost_type));
  }
  
  /// loop on each quadrature point
  for (; traction_it != traction_end;
       ++traction_it, ++opening_it, ++opening_prec_it, 
	 ++contact_traction_it, ++contact_opening_it, ++normal_it,
	 ++sigma_c_it, ++delta_n_max_it, ++delta_t_max_it,
	 ++delta_c_it, ++damage_n_it, ++damage_t_it,
	 ++insertion_stress_it, ++reduction_penalty_it) {

    Real normal_opening_norm, tangential_opening_norm;
    bool penetration;
    Real current_penalty = 0.;
    Real delta_c2_R2 = *delta_c_it * (*delta_c_it) / (R * R);
    
    /// compute normal and tangential opening vectors
    normal_opening_norm = opening_it->dot(*normal_it);
    Vector<Real> normal_opening  = *normal_it;
    normal_opening *= normal_opening_norm;

    Vector<Real> tangential_opening  = *opening_it;
    tangential_opening -= normal_opening;
    tangential_opening_norm = tangential_opening.norm();
  
    /// compute effective opening displacement
    Real delta_n = tangential_opening_norm * tangential_opening_norm * this->beta2_kappa2;
    Real delta_t = tangential_opening_norm * tangential_opening_norm * this->beta2_kappa2;

    penetration = normal_opening_norm < 0.0;
    //    if (this->contact_after_breaking == false && Math::are_float_equal(damage, 1.))
    //  penetration = false;

    /** 
     * If during the convergence loop a cohesive element continues to
     * jumps from penetration to opening, and convergence is not
     * reached, its penalty parameter will be reduced in the
     * recomputation of the same incremental step. Recompute is set
     * equal to true when convergence is not reached in the
     * solveStepCohesive function and the execution of the program
     * goes back to the main file where the variable load_reduction
     * is set equal to true.
     */
    Real normal_opening_prec_norm = opening_prec_it->dot(*normal_it);
    Vector<Real> normal_opening_prec  = *normal_it;
    normal_opening_prec *= normal_opening_prec_norm;
    if (!this->model->isExplicit())  // && !this->recompute)
      if ((normal_opening_prec_norm * normal_opening_norm) < 0.0)
	*reduction_penalty_it = true;

    *opening_prec_it = *opening_it;
    
    if (penetration) {
      if (this->recompute && *reduction_penalty_it){
	/// the penalty parameter is locally reduced
	current_penalty = this->penalty / 10.;
      }
      else
	current_penalty = this->penalty;

      /// use penalty coefficient in case of penetration
      *contact_traction_it = normal_opening;
      *contact_traction_it *= current_penalty;
      *contact_opening_it = normal_opening;

      /// don't consider penetration contribution for delta
      *opening_it = tangential_opening;
      normal_opening.clear();

    }
    else {
      delta_n += normal_opening_norm * normal_opening_norm;
      delta_t += normal_opening_norm * normal_opening_norm * delta_c2_R2;
      contact_traction_it->clear();
      contact_opening_it->clear();
    }

    delta_n = std::sqrt(delta_n);
    delta_t = std::sqrt(delta_t);

    /// update maximum displacement and damage
    *delta_n_max_it = std::max(*delta_n_max_it, delta_n);
    *damage_n_it = std::min(*delta_n_max_it / *delta_c_it, Real(1.));

    *delta_t_max_it = std::max(*delta_t_max_it, delta_t);
    *damage_t_it = std::min(*delta_t_max_it / *delta_c_it, Real(1.));

    Vector<Real> normal_traction(spatial_dimension);
    Vector<Real> shear_traction(spatial_dimension);

    /// NORMAL TRACTIONS
    if (Math::are_float_equal(*damage_n_it, 1.))
      normal_traction.clear();
    else if (Math::are_float_equal(*damage_n_it, 0.)) {
      if (penetration)
	normal_traction.clear();
      else
	normal_traction = *insertion_stress_it;
    }
    else {
      // the following formulation holds both in loading and in
      // unloading-reloading
      normal_traction = normal_opening;

      AKANTU_DEBUG_ASSERT(*delta_n_max_it != 0.,
			  "Division by zero, tolerance might be too low");

      normal_traction *= *sigma_c_it / (*delta_n_max_it) * (1. - *damage_n_it);
    }

    /// SHEAR TRACTIONS
    if (Math::are_float_equal(*damage_t_it, 1.))
      shear_traction.clear();
    else if (Math::are_float_equal(*damage_t_it, 0.)) {
      shear_traction.clear();
    }
    else {
      shear_traction = tangential_opening;

      AKANTU_DEBUG_ASSERT(*delta_t_max_it != 0.,
			  "Division by zero, tolerance might be too low");

      shear_traction *= *sigma_c_it * this->beta2_kappa / (*delta_t_max_it) * (1. - *damage_t_it);
    }

    *traction_it = normal_traction;
    *traction_it += shear_traction;

  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveLinearUncoupled<spatial_dimension>::checkDeltaMax(GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  
  /**
  * This function set a predefined value to the parameter
  * delta_max_prev of the elements that have been inserted in the
  * last loading step for which convergence has not been
  * reached. This is done before reducing the loading and re-doing
  * the step.  Otherwise, the updating of delta_max_prev would be
  * done with reference to the non-convergent solution. In this
  * function also other variables related to the contact and
  * friction behavior are correctly set.
  */

  Mesh & mesh = this->fem_cohesive->getMesh();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension,
					  ghost_type, _ek_cohesive);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension,
						ghost_type, _ek_cohesive);

  /**
  * the variable "recompute" is set to true to activate the
  * procedure that reduces the penalty parameter for
  * compression. This procedure is available only during the phase of
  * load_reduction, that has to be set in the main file. The
  * penalty parameter will be reduced only for the elements having
  * reduction_penalty = true.
  */
  this->recompute = true;

  for(; it != last_type; ++it) {
    Array<UInt> & elem_filter = this->element_filter(*it, ghost_type);

    UInt nb_element = elem_filter.getSize();
    if (nb_element == 0) continue;

    ElementType el_type = *it;

    //    std::cout << "element type: " << el_type << std::endl;

    /// define iterators
    Array<Real>::scalar_iterator delta_n_max_it =
      delta_n_max(el_type, ghost_type).begin();

    Array<Real>::scalar_iterator delta_n_max_end =
      delta_n_max(el_type, ghost_type).end();

    Array<Real>::scalar_iterator delta_n_max_prev_it =
      delta_n_max.previous(el_type, ghost_type).begin();

    Array<Real>::scalar_iterator delta_t_max_it =
      delta_t_max(el_type, ghost_type).begin();

    Array<Real>::scalar_iterator delta_t_max_prev_it =
      delta_t_max.previous(el_type, ghost_type).begin();

    Array<Real>::scalar_iterator delta_c_it =
     this-> delta_c_eff(el_type, ghost_type).begin();

    Array<Real>::vector_iterator opening_prec_it =
      this->opening_prec(el_type, ghost_type).begin(spatial_dimension);

    Array<Real>::vector_iterator opening_prec_prev_it =
      this->opening_prec.previous(el_type, ghost_type).begin(spatial_dimension);

  Int it = 0;

    /// loop on each quadrature point
    for (; delta_n_max_it != delta_n_max_end;
         ++delta_n_max_it, ++delta_t_max_it, ++delta_c_it,
	   ++delta_n_max_prev_it, ++delta_t_max_prev_it, 
	   ++opening_prec_it, ++opening_prec_prev_it) {

      ++it;

      if (*delta_n_max_prev_it == 0){
        /// elements inserted in the last incremental step, that did
        /// not converge
        *delta_n_max_it = *delta_c_it / 1000.;
      }
      else
        /// elements introduced in previous incremental steps, for
        /// which a correct value of delta_max_prev already exists
        *delta_n_max_it = *delta_n_max_prev_it;

      if (*delta_t_max_prev_it == 0){
        *delta_t_max_it = *delta_c_it * this->kappa / this->beta / 1000.;
      }
      else
        *delta_t_max_it = *delta_t_max_prev_it;

      /// in case convergence is not reached, set opening_prec to the
      /// value referred to the previous incremental step
      *opening_prec_it = *opening_prec_prev_it;
    }
  }
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveLinearUncoupled<spatial_dimension>::computeTangentTraction(const ElementType & el_type,
                                                                       Array<Real> & tangent_matrix,
                                                                       const Array<Real> & normal,
                                                                       GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  /// define iterators
  Array<Real>::matrix_iterator tangent_it =
    tangent_matrix.begin(spatial_dimension, spatial_dimension);

  Array<Real>::matrix_iterator tangent_end =
    tangent_matrix.end(spatial_dimension, spatial_dimension);

  Array<Real>::const_vector_iterator normal_it =
    normal.begin(spatial_dimension);

  Array<Real>::vector_iterator opening_it =
    this->opening(el_type, ghost_type).begin(spatial_dimension);

  /// NB: delta_max_it points on delta_max_previous, i.e. the
  /// delta_max related to the solution of the previous incremental
  /// step
  Array<Real>::scalar_iterator delta_n_max_it =
   this->delta_n_max.previous(el_type, ghost_type).begin();

  Array<Real>::scalar_iterator delta_t_max_it =
    this->delta_t_max.previous(el_type, ghost_type).begin();

  Array<Real>::scalar_iterator sigma_c_it =
    this->sigma_c_eff(el_type, ghost_type).begin();

  Array<Real>::scalar_iterator delta_c_it =
    this->delta_c_eff(el_type, ghost_type).begin();

  Array<Real>::vector_iterator contact_opening_it =
    this->contact_opening(el_type, ghost_type).begin(spatial_dimension);

  Array<bool>::scalar_iterator reduction_penalty_it =
    this->reduction_penalty(el_type, ghost_type).begin();

  Vector<Real> normal_opening(spatial_dimension);
  Vector<Real> tangential_opening(spatial_dimension);

  for (; tangent_it != tangent_end;
       ++tangent_it, ++normal_it, ++opening_it, ++sigma_c_it,
	 ++delta_c_it, ++delta_n_max_it, ++delta_t_max_it,
	 ++contact_opening_it, ++reduction_penalty_it) {

    Real normal_opening_norm, tangential_opening_norm;
    bool penetration;
    Real current_penalty = 0.;
    Real delta_c2_R2 = *delta_c_it * (*delta_c_it) / (R * R);

    /**
     * During the update of the residual the interpenetrations are
     * stored in the array "contact_opening", therefore, in the case
     * of penetration, in the array "opening" there are only the
     * tangential components.
     */
    *opening_it += *contact_opening_it;

    /// compute normal and tangential opening vectors
    normal_opening_norm = opening_it->dot(*normal_it);
    Vector<Real> normal_opening = *normal_it;
    normal_opening *= normal_opening_norm;

    Vector<Real> tangential_opening = *opening_it;
    tangential_opening -= normal_opening;
    tangential_opening_norm = tangential_opening.norm();

    Real delta_n = tangential_opening_norm * tangential_opening_norm * this->beta2_kappa2;
    Real delta_t = tangential_opening_norm * tangential_opening_norm * this->beta2_kappa2;

    penetration = normal_opening_norm < 0.0;
    //    if (this->contact_after_breaking == false && Math::are_float_equal(damage, 1.))
    //  penetration = false;

    Real derivative = 0; // derivative = d(t/delta)/ddelta
    Real T = 0;
    
    /// TANGENT STIFFNESS FOR NORMAL TRACTIONS
    Matrix<Real> n_outer_n(spatial_dimension, spatial_dimension);
    n_outer_n.outerProduct(*normal_it, *normal_it);

    if (penetration){
      if (this->recompute && *reduction_penalty_it)
	current_penalty = this->penalty / 10.;
      else
	current_penalty = this->penalty;

      /// stiffness in compression given by the penalty parameter
      *tangent_it = n_outer_n;
      *tangent_it *= current_penalty;
  
      *opening_it = tangential_opening;
      normal_opening.clear();
    }
    else{
      delta_n += normal_opening_norm * normal_opening_norm;
      delta_n = std::sqrt(delta_n);

      delta_t += normal_opening_norm * normal_opening_norm * delta_c2_R2;
      
      /**
       * Delta has to be different from 0 to have finite values of
       * tangential stiffness.  At the element insertion, delta =
       * 0. Therefore, a fictictious value is defined, for the
       * evaluation of the first value of K.
       */
      if (delta_n < Math::getTolerance())
	delta_n = *delta_c_it / 1000.;

      // loading
      if (delta_n >= *delta_n_max_it){
	if (delta_n <= *delta_c_it){
	  derivative = - (*sigma_c_it) / (delta_n * delta_n);
	  T = *sigma_c_it * (1 - delta_n / *delta_c_it);
	} else {
	  derivative = 0.;
	  T = 0.;
	}
	// unloading-reloading
      } else if (delta_n < *delta_n_max_it){
	Real T_max = *sigma_c_it * (1 - *delta_n_max_it / *delta_c_it);
	derivative = 0.;
	T = T_max / *delta_n_max_it * delta_n;
      }

      /// computation of the derivative of the constitutive law (dT/ddelta)
      Matrix<Real> nn(n_outer_n);
      nn *= T / delta_n;

      Vector<Real> Delta_tilde(normal_opening);
      Delta_tilde *= (1. - this->beta2_kappa2);
      Vector<Real> mm(*opening_it);
      mm *= this->beta2_kappa2;
      Delta_tilde += mm;

      Vector<Real> Delta_hat(normal_opening);
      Matrix<Real> prov(spatial_dimension, spatial_dimension);
      prov.outerProduct(Delta_hat, Delta_tilde);
      prov *= derivative / delta_n;
      prov += nn;

      *tangent_it = prov.transpose();
    }

    derivative = 0.;
    T = 0.;

    /// TANGENT STIFFNESS FOR SHEAR TRACTIONS
    delta_t = std::sqrt(delta_t);

    /**
     * Delta has to be different from 0 to have finite values of
     * tangential stiffness.  At the element insertion, delta =
     * 0. Therefore, a fictictious value is defined, for the
     * evaluation of the first value of K.
     */
    if (delta_t < Math::getTolerance())
      delta_t = *delta_c_it / 1000.;

    // loading
    if (delta_t >= *delta_t_max_it){
      if (delta_t <= *delta_c_it){
	derivative = - (*sigma_c_it) / (delta_t * delta_t);
	T = *sigma_c_it * (1 - delta_t / *delta_c_it);
      } else {
	derivative = 0.;
	T = 0.;
      }
      // unloading-reloading
    } else if (delta_t < *delta_t_max_it){
      Real T_max = (*sigma_c_it) * (1 - *delta_t_max_it / *delta_c_it);
      derivative = 0.;
      T = T_max / *delta_t_max_it * delta_t;
    }

    /// computation of the derivative of the constitutive law (dT/ddelta)
    Matrix<Real> I(spatial_dimension, spatial_dimension);
    I.eye();
    Matrix<Real> nn(n_outer_n);
    I -= nn;
    I *= T / delta_t;
    
    Vector<Real> Delta_tilde(normal_opening);
    Delta_tilde *= (delta_c2_R2 - this->beta2_kappa2);
    Vector<Real> mm(*opening_it);
    mm *= this->beta2_kappa2;
    Delta_tilde += mm;

    Vector<Real> Delta_hat(tangential_opening);
    Delta_hat *= this->beta2_kappa;
    Matrix<Real> prov(spatial_dimension, spatial_dimension);
    prov.outerProduct(Delta_hat, Delta_tilde);
    prov *= derivative / delta_t;
    prov += I;

    //    Matrix<Real> prov_t = prov.transpose();
    *tangent_it += prov.transpose();
  }
    
  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */


INSTANTIATE_MATERIAL(MaterialCohesiveLinearUncoupled);

__END_AKANTU__
