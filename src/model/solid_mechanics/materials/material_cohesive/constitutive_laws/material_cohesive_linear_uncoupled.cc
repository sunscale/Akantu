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
  delta_nt_max("delta_nt_max", *this),
  delta_st_max("delta_st_max", *this),
  damage_nt("damage_nt", *this),
  damage_st("damage_st", *this){
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

  delta_nt_max.initialize(1);
  delta_st_max.initialize(1);
  damage_nt.initialize(1);
  damage_st.initialize(1);
  
  delta_nt_max.initializeHistory();
  delta_st_max.initializeHistory();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveLinearUncoupled<spatial_dimension>::computeTraction(const Array<Real> & normal,
                                                                ElementType el_type,
                                                                GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  delta_nt_max.resize();
  delta_st_max.resize();
  damage_nt.resize();
  damage_st.resize();
  
  /// define iterators
  Array<Real>::vector_iterator traction_it =
    this->tractions(el_type, ghost_type).begin(spatial_dimension);

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

  Array<Real>::vector_iterator traction_end =
    this->tractions(el_type, ghost_type).end(spatial_dimension);

  Array<Real>::scalar_iterator sigma_c_it =
    this->sigma_c_eff(el_type, ghost_type).begin();

  Array<Real>::scalar_iterator delta_nt_max_it =
    this->delta_nt_max(el_type, ghost_type).begin();

  Array<Real>::scalar_iterator delta_st_max_it =
    this->delta_st_max(el_type, ghost_type).begin();

  Array<Real>::scalar_iterator delta_nt_max_prev_it =
    this->delta_nt_max.previous(el_type, ghost_type).begin();

  Array<Real>::scalar_iterator delta_st_max_prev_it =
    this->delta_st_max.previous(el_type, ghost_type).begin();

  Array<Real>::scalar_iterator delta_c_it =
    this->delta_c_eff(el_type, ghost_type).begin();

  Array<Real>::scalar_iterator damage_nt_it =
    this->damage_nt(el_type, ghost_type).begin();

  Array<Real>::scalar_iterator damage_st_it =
    this->damage_st(el_type, ghost_type).begin();

  Array<Real>::vector_iterator insertion_stress_it =
    this->insertion_stress(el_type, ghost_type).begin(spatial_dimension);

  Array<bool>::scalar_iterator reduction_penalty_it =
    this->reduction_penalty(el_type, ghost_type).begin();

  Vector<Real> normal_opening(spatial_dimension);
  Vector<Real> tangential_opening(spatial_dimension);

  if (! this->model->isExplicit()){
    this->delta_nt_max(el_type, ghost_type).copy(this->delta_nt_max.previous(el_type, ghost_type));
    this->delta_st_max(el_type, ghost_type).copy(this->delta_st_max.previous(el_type, ghost_type));
  }
  
  /// loop on each quadrature point
  for (; traction_it != traction_end;
       ++traction_it, ++opening_it, ++opening_prec_it,
	 ++normal_it, ++sigma_c_it, ++delta_c_it, ++delta_nt_max_it, 
	 ++delta_st_max_it, ++damage_nt_it, ++damage_st_it,
	 ++contact_traction_it, ++contact_opening_it, ++insertion_stress_it,
	 ++reduction_penalty_it) {

    Real normal_opening_norm, tangential_opening_norm;
    bool penetration;
    Real current_penalty = 0.;
    Real delta_c_s = this->kappa / this->beta * (*delta_c_it);
    Real delta_c_s2_R2 = delta_c_s * delta_c_s / (R * R);
    
    /// compute normal and tangential opening vectors
    normal_opening_norm = opening_it->dot(*normal_it);
    Vector<Real> normal_opening  = *normal_it;
    normal_opening *= normal_opening_norm;

    Vector<Real> tangential_opening  = *opening_it;
    tangential_opening -= normal_opening;
    tangential_opening_norm = tangential_opening.norm();
  
    /// compute effective opening displacement
    Real delta_nt = tangential_opening_norm * tangential_opening_norm * this->beta2_kappa2;
    Real delta_st = tangential_opening_norm * tangential_opening_norm;

    penetration = normal_opening_norm < 0.0;
    //    if (this->contact_after_breaking == false && Math::are_float_equal(damage, 1.))
    //  penetration = false;

    /** 
     * if during the convergence loop a cohesive element continues to
     * jumps from penetration to opening, and convergence is not
     * reached, its penalty parameter will be reduced in the
     * recomputation of the same incremental step. recompute is set
     * equal to true when convergence is not reached in the
     * solveStepCohesive function and the execution of the program
     * goes back to the main file where the variable load_reduction
     * is set equal to true.
     */
    Real normal_opening_prec_norm = opening_prec_it->dot(*normal_it);
    *opening_prec_it = *opening_it;

    if (!this->model->isExplicit() && !this->recompute)
      if ((normal_opening_prec_norm * normal_opening_norm) < 0.0) {
	*reduction_penalty_it = true;
      }

    if (penetration) {
      if (this->recompute && *reduction_penalty_it){
	/// the penalty parameter is locally reduced
	current_penalty = this->penalty / 100.;
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
      delta_nt += normal_opening_norm * normal_opening_norm;
      delta_st += normal_opening_norm * normal_opening_norm * delta_c_s2_R2;
      contact_traction_it->clear();
      contact_opening_it->clear();
    }

    delta_nt = std::sqrt(delta_nt);
    delta_st = std::sqrt(delta_st);

    /// update maximum displacement and damage
    *delta_nt_max_it = std::max(*delta_nt_max_it, delta_nt);
    *damage_nt_it = std::min(*delta_nt_max_it / *delta_c_it, Real(1.));

    *delta_st_max_it = std::max(*delta_st_max_it, delta_st);
    *damage_st_it = std::min(*delta_st_max_it / delta_c_s, Real(1.));

    /**
     * Compute traction @f$ \mathbf{T} = \left(
     * \frac{\beta^2}{\kappa} \Delta_t \mathbf{t} + \Delta_n
     * \mathbf{n} \right) \frac{\sigma_c}{\delta} \left( 1-
     * \frac{\delta}{\delta_c} \right)@f$
     */

    Vector<Real> traction_n(spatial_dimension);
    Vector<Real> traction_s(spatial_dimension);

    /// NORMAL TRACTIONS
    if (Math::are_float_equal(*damage_nt_it, 1.))
      traction_n.clear();
    else if (Math::are_float_equal(*damage_nt_it, 0.)) {
      if (penetration)
	traction_n.clear();
      else
	traction_n = *insertion_stress_it;
    }
    else {
      traction_n = normal_opening;

      AKANTU_DEBUG_ASSERT(*delta_nt_max_it != 0.,
			  "Division by zero, tolerance might be too low");

      traction_n *= *sigma_c_it / (*delta_nt_max_it) * (1. - *damage_nt_it);
    }

    /// SHEAR TRACTIONS
    if (Math::are_float_equal(*damage_st_it, 1.))
      traction_s.clear();
    else if (Math::are_float_equal(*damage_st_it, 0.)) {
      traction_s.clear();
    }
    else {
      traction_s = tangential_opening;

      AKANTU_DEBUG_ASSERT(*delta_st_max_it != 0.,
			  "Division by zero, tolerance might be too low");

      traction_s *= *sigma_c_it * this->beta / (*delta_st_max_it) * (1. - *damage_st_it);
    }

    *traction_it = traction_n + traction_s;

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
  * procedure that reduce the penalty parameter for
  * compression. This procedure is set true only during the phase of
  * load_reduction, that has to be set in the maiin file. The
  * penalty parameter will be reduced only for the elements having
  * reduction_penalty = true.
  */
  this->recompute = true;

  for(; it != last_type; ++it) {
    Array<UInt> & elem_filter = this->element_filter(*it, ghost_type);

    UInt nb_element = elem_filter.getSize();
    if (nb_element == 0) continue;

    ElementType el_type = *it;

    /// define iterators
    Array<Real>::scalar_iterator delta_nt_max_it =
      delta_nt_max(el_type, ghost_type).begin();

    Array<Real>::scalar_iterator delta_nt_max_end =
      delta_nt_max(el_type, ghost_type).end();

    Array<Real>::scalar_iterator delta_nt_max_prev_it =
      delta_nt_max.previous(el_type, ghost_type).begin();

    Array<Real>::scalar_iterator delta_st_max_it =
      delta_st_max(el_type, ghost_type).begin();

    Array<Real>::scalar_iterator delta_st_max_prev_it =
      delta_st_max.previous(el_type, ghost_type).begin();

    Array<Real>::scalar_iterator delta_c_it =
     this-> delta_c_eff(el_type, ghost_type).begin();

    Array<Real>::vector_iterator opening_prec_it =
      this->opening_prec(el_type, ghost_type).begin(spatial_dimension);

    Array<Real>::vector_iterator opening_prec_prev_it =
      this->opening_prec.previous(el_type, ghost_type).begin(spatial_dimension);

    /// loop on each quadrature point
    for (; delta_nt_max_it != delta_nt_max_end;
         ++delta_nt_max_it, ++delta_st_max_it, ++delta_c_it,
	   ++delta_nt_max_prev_it, ++delta_st_max_prev_it, 
	   ++opening_prec_it, ++opening_prec_prev_it) {

      if (*delta_nt_max_prev_it == 0)
        /// elements inserted in the last incremental step, that did
        /// not converge
        *delta_nt_max_it = *delta_c_it / 1000;
      else
        /// elements introduced in previous incremental steps, for
        /// which a correct value of delta_max_prev already exists
        *delta_nt_max_it = *delta_nt_max_prev_it;

      if (*delta_st_max_prev_it == 0)
        *delta_st_max_it = *delta_c_it * this->kappa / this->beta / 1000;
      else
        *delta_st_max_it = *delta_st_max_prev_it;

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
  Array<Real>::scalar_iterator delta_nt_max_it =
   this->delta_nt_max.previous(el_type, ghost_type).begin();

  Array<Real>::scalar_iterator delta_st_max_it =
    this->delta_st_max.previous(el_type, ghost_type).begin();

  Array<Real>::scalar_iterator sigma_c_it =
    this->sigma_c_eff(el_type, ghost_type).begin();

  Array<Real>::scalar_iterator delta_c_it =
    this->delta_c_eff(el_type, ghost_type).begin();

  //  Array<Real>::scalar_iterator damage_it =
  //    damage(el_type, ghost_type).begin();

  Array<Real>::vector_iterator contact_opening_it =
    this->contact_opening(el_type, ghost_type).begin(spatial_dimension);

  Array<bool>::scalar_iterator reduction_penalty_it =
    this->reduction_penalty(el_type, ghost_type).begin();

  Vector<Real> normal_opening(spatial_dimension);
  Vector<Real> tangential_opening(spatial_dimension);

  for (; tangent_it != tangent_end;
       ++tangent_it, ++normal_it, ++opening_it, ++sigma_c_it,
	 ++delta_c_it, ++delta_nt_max_it, ++delta_st_max_it,
	 ++contact_opening_it, ++reduction_penalty_it) {

    Real normal_opening_norm, tangential_opening_norm;
    bool penetration;
    Real current_penalty = 0.;
    Real delta_c_s = this->kappa / this->beta * (*delta_c_it);
    Real delta_c_s2_R2 = delta_c_s * delta_c_s / (R * R);

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

    Real delta_nt = tangential_opening_norm * tangential_opening_norm * this->beta2_kappa2;
    Real delta_st = tangential_opening_norm * tangential_opening_norm;

    penetration = normal_opening_norm < 0.0;
    //    if (this->contact_after_breaking == false && Math::are_float_equal(damage, 1.))
    //  penetration = false;

    Real derivative = 0; // derivative = d(t/delta)/ddelta
    Real t = 0;
    //    Matrix<Real> tangent_n(spatial_dimension, spatial_dimension);
    //    Matrix<Real> tangent_s(spatial_dimension, spatial_dimension);
    
    /// TANGENT STIFFNESS FOR NORMAL TRACTIONS
    Matrix<Real> n_outer_n(spatial_dimension, spatial_dimension);
    n_outer_n.outerProduct(*normal_it, *normal_it);

    if (penetration){
      if (this->recompute && *reduction_penalty_it)
	current_penalty = this->penalty / 100.;
      else
	current_penalty = this->penalty;

      /// stiffness in compression given by the penalty parameter
      *tangent_it = n_outer_n;
      *tangent_it *= current_penalty;
  
      *opening_it = tangential_opening;
      normal_opening.clear();
    }
    else{
      delta_nt += normal_opening_norm * normal_opening_norm;
      delta_nt = std::sqrt(delta_nt);

      delta_st += normal_opening_norm * normal_opening_norm * delta_c_s2_R2;
      
      /**
       * Delta has to be different from 0 to have finite values of
       * tangential stiffness.  At the element insertion, delta =
       * 0. Therefore, a fictictious value is defined, for the
       * evaluation of the first value of K.
       */
      if (delta_nt < Math::getTolerance())
	delta_nt = *delta_c_it / 1000.;

      if (delta_nt >= *delta_nt_max_it){
	if (delta_nt <= *delta_c_it){
	  derivative = - (*sigma_c_it) / (delta_nt * delta_nt);
	  t = *sigma_c_it * (1 - delta_nt / *delta_c_it);
	} else {
	  derivative = 0.;
	  t = 0.;
	}
      } else if (delta_nt < *delta_nt_max_it){
	Real tmax = *sigma_c_it * (1 - *delta_nt_max_it / *delta_c_it);
	t = tmax / *delta_nt_max_it * delta_nt;
      }

      /// computation of the derivative of the constitutive law (dT/ddelta)
      Matrix<Real> nn(n_outer_n);
      nn *= t / delta_nt;

      Vector<Real> t_tilde(normal_opening);
      t_tilde *= (1. - this->beta2_kappa2);
      Vector<Real> mm(*opening_it);
      mm *= this->beta2_kappa2;
      t_tilde += mm;

      Vector<Real> t_hat(normal_opening);
      Matrix<Real> prov(spatial_dimension, spatial_dimension);
      prov.outerProduct(t_hat, t_tilde);
      prov *= derivative / delta_nt;
      prov += nn;

      *tangent_it = prov.transpose();
    }

    derivative = 0.;
    t = 0.;
    /// TANGENT STIFFNESS FOR SHEAR TRACTIONS
    delta_st = std::sqrt(delta_st);

    /**
     * Delta has to be different from 0 to have finite values of
     * tangential stiffness.  At the element insertion, delta =
     * 0. Therefore, a fictictious value is defined, for the
     * evaluation of the first value of K.
     */
    if (delta_st < Math::getTolerance())
      delta_st = delta_c_s / 1000.;

    if (delta_st >= *delta_st_max_it){
      if (delta_st <= delta_c_s){
	derivative = - this->beta * (*sigma_c_it) / (delta_st * delta_st);
	t = this->beta * (*sigma_c_it) * (1 - delta_st / delta_c_s);
      } else {
	derivative = 0.;
	t = 0.;
      }
    } else if (delta_st < *delta_st_max_it){
      Real tmax = this->beta * (*sigma_c_it) * (1 - *delta_st_max_it / delta_c_s);
      t = tmax / *delta_st_max_it * delta_st;
    }

    /// computation of the derivative of the constitutive law (dT/ddelta)
    Matrix<Real> I(spatial_dimension, spatial_dimension);
    I.eye();
    Matrix<Real> nn(n_outer_n);
    I -= nn;
    I *= t / delta_st;
    
    Vector<Real> t_tilde(normal_opening);
    t_tilde *= (delta_c_s2_R2 - 1.);
    Vector<Real> mm(*opening_it);
    t_tilde += mm;

    Vector<Real> t_hat(tangential_opening);
    Matrix<Real> prov(spatial_dimension, spatial_dimension);
    prov.outerProduct(t_hat, t_tilde);
    prov *= derivative / delta_st;
    prov += I;

    Matrix<Real> prov_t = prov.transpose();
    *tangent_it += prov_t;
  }
    
  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */


INSTANTIATE_MATERIAL(MaterialCohesiveLinearUncoupled);

__END_AKANTU__
