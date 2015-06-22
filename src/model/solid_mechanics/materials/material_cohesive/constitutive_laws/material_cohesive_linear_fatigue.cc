/**
 * @file   material_cohesive_linear_fatigue.cc
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Thu Feb 19 14:40:57 2015
 *
 * @brief  See material_cohesive_linear_fatigue.hh for information
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
#include "material_cohesive_linear_fatigue.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
MaterialCohesiveLinearFatigue<spatial_dimension>
::MaterialCohesiveLinearFatigue(SolidMechanicsModel & model,
				const ID & id) :
  MaterialCohesiveLinear<spatial_dimension>(model, id),
  delta_prec("delta_prec", *this),
  K_plus("K_plus", *this),
  K_minus("K_minus", *this),
  T_1d("T_1d", *this),
  switches("switches", *this),
  delta_dot_prec("delta_dot_prec", *this) {

  this->registerParam("delta_f", delta_f, -1. ,
		      _pat_parsable | _pat_readable,
		      "delta_f");

  this->registerParam("progressive_delta_f", progressive_delta_f, false,
		      _pat_parsable | _pat_readable,
		      "Whether or not delta_f is equal to delta_max");

  this->registerParam("count_switches", count_switches, false,
		      _pat_parsable | _pat_readable,
		      "Count the opening/closing switches per element");
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialCohesiveLinearFatigue<spatial_dimension>::initMaterial() {
  MaterialCohesiveLinear<spatial_dimension>::initMaterial();

  // check that delta_f has a proper value or assign a defaul value
  if (delta_f < 0) delta_f = this->delta_c_eff;
  else if (delta_f < this->delta_c_eff)
    AKANTU_DEBUG_ERROR("Delta_f must be greater or equal to delta_c");

  delta_prec.initialize(1);
  K_plus.initialize(1);
  K_minus.initialize(1);
  T_1d.initialize(1);

  if (count_switches) {
    switches.initialize(1);
    delta_dot_prec.initialize(1);
  }
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveLinearFatigue<spatial_dimension>
::computeTraction(const Array<Real> & normal,
		  ElementType el_type,
		  GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  /// define iterators
  Array<Real>::vector_iterator traction_it =
    this->tractions(el_type, ghost_type).begin(spatial_dimension);

  Array<Real>::vector_iterator opening_it =
    this->opening(el_type, ghost_type).begin(spatial_dimension);

  Array<Real>::vector_iterator contact_traction_it =
    this->contact_tractions(el_type, ghost_type).begin(spatial_dimension);

  Array<Real>::vector_iterator contact_opening_it =
    this->contact_opening(el_type, ghost_type).begin(spatial_dimension);

  Array<Real>::const_vector_iterator normal_it = normal.begin(spatial_dimension);

  Array<Real>::vector_iterator traction_end =
    this->tractions(el_type, ghost_type).end(spatial_dimension);

  Array<Real>::scalar_iterator sigma_c_it =
    this->sigma_c_eff(el_type, ghost_type).begin();

  Array<Real>::scalar_iterator delta_max_it =
    this->delta_max(el_type, ghost_type).begin();

  Array<Real>::scalar_iterator delta_c_it =
    this->delta_c_eff(el_type, ghost_type).begin();

  Array<Real>::scalar_iterator damage_it = this->damage(el_type, ghost_type).begin();

  Array<Real>::vector_iterator insertion_stress_it =
    this->insertion_stress(el_type, ghost_type).begin(spatial_dimension);

  Array<Real>::scalar_iterator delta_prec_it = delta_prec(el_type, ghost_type).begin();
  Array<Real>::scalar_iterator K_plus_it = K_plus(el_type, ghost_type).begin();
  Array<Real>::scalar_iterator K_minus_it = K_minus(el_type, ghost_type).begin();
  Array<Real>::scalar_iterator T_1d_it = T_1d(el_type, ghost_type).begin();

  Array<UInt>::scalar_iterator switches_it;
  Array<Real>::scalar_iterator delta_dot_prec_it;

  if (count_switches) {
    switches_it = switches(el_type, ghost_type).begin();
    delta_dot_prec_it = delta_dot_prec(el_type, ghost_type).begin();
  }

  Real * memory_space = new Real[2*spatial_dimension];
  Vector<Real> normal_opening(memory_space, spatial_dimension);
  Vector<Real> tangential_opening(memory_space + spatial_dimension,
				  spatial_dimension);

  Real tolerance = Math::getTolerance();

  /// loop on each quadrature point
  for (; traction_it != traction_end;
       ++traction_it, ++opening_it, ++normal_it, ++sigma_c_it,
	 ++delta_max_it, ++delta_c_it, ++damage_it, ++contact_traction_it,
	 ++insertion_stress_it, ++contact_opening_it, ++delta_prec_it,
	 ++K_plus_it, ++K_minus_it, ++T_1d_it) {

    /// compute normal and tangential opening vectors
    Real normal_opening_norm = opening_it->dot(*normal_it);
    normal_opening  = (*normal_it);
    normal_opening *= normal_opening_norm;

    tangential_opening  = *opening_it;
    tangential_opening -=  normal_opening;

    Real tangential_opening_norm = tangential_opening.norm();

    /**
     * compute effective opening displacement
     * @f$ \delta = \sqrt{
     * \frac{\beta^2}{\kappa^2} \Delta_t^2 + \Delta_n^2 } @f$
     */
    Real delta = tangential_opening_norm * tangential_opening_norm * this->beta2_kappa2;

    bool penetration = normal_opening_norm < -tolerance;
    if (this->contact_after_breaking == false && Math::are_float_equal(*damage_it, 1.))
      penetration = false;

    if (penetration) {
      /// use penalty coefficient in case of penetration
      *contact_traction_it = normal_opening;
      *contact_traction_it *= this->penalty;
      *contact_opening_it = normal_opening;
      /// don't consider penetration contribution for delta
      *opening_it = tangential_opening;
      normal_opening.clear();
    }
    else {
      delta += normal_opening_norm * normal_opening_norm;
      contact_traction_it->clear();
      contact_opening_it->clear();
    }

    delta = std::sqrt(delta);

    /**
     * Compute traction @f$ \mathbf{T} = \left(
     * \frac{\beta^2}{\kappa} \Delta_t \mathbf{t} + \Delta_n
     * \mathbf{n} \right) \frac{\sigma_c}{\delta} \left( 1-
     * \frac{\delta}{\delta_c} \right)@f$
     */


    // update maximum displacement and damage
    *delta_max_it = std::max(delta, *delta_max_it);
    *damage_it = std::min(*delta_max_it / *delta_c_it, 1.);

    Real delta_dot = delta - *delta_prec_it;

    // count switches if asked
    if (count_switches) {
      if ( (delta_dot > 0. && *delta_dot_prec_it <= 0.) ||
	   (delta_dot < 0. && *delta_dot_prec_it >= 0.) )
	++(*switches_it);

      *delta_dot_prec_it = delta_dot;
      ++delta_dot_prec_it;
      ++switches_it;
    }

    // set delta_f equal to delta_max if desired
    if (progressive_delta_f) delta_f = *delta_max_it;

    // broken element case
    if (Math::are_float_equal(*damage_it, 1.))
      traction_it->clear();
    // just inserted element case
    else if (Math::are_float_equal(*damage_it, 0.)) {
      if (penetration)
	traction_it->clear();
      else
	*traction_it = *insertion_stress_it;
      // initialize the 1d traction to sigma_c
      *T_1d_it = *sigma_c_it;
    }
    // normal case
    else {
      // if element is closed then there are zero tractions
      if (delta <= tolerance)
	traction_it->clear();
      // otherwise compute new tractions if the new delta is different
      // than the previous one
      else if (std::abs(delta_dot) > tolerance) {
	// loading case
	if (delta_dot > 0.) {
	  // equation (4) of the article
	  *K_plus_it *= 1. - delta_dot / delta_f;
	  // equivalent to equation (2) of the article
	  *T_1d_it += *K_plus_it * delta_dot;

	  // in case of reloading, traction must not exceed that of the
	  // envelop of the cohesive law
	  Real max_traction = *sigma_c_it * (1 - delta / *delta_c_it);
	  bool max_traction_exceeded = *T_1d_it > max_traction;
	  if (max_traction_exceeded) *T_1d_it = max_traction;

	  // equation (3) of the article
	  *K_minus_it = *T_1d_it / delta;

	  // if the traction is following the cohesive envelop, then
	  // K_plus has to be reset
	  if (max_traction_exceeded) *K_plus_it = *K_minus_it;
	}
	// unloading case
	else {
	  // equation (4) of the article
	  *K_plus_it += (*K_plus_it - *K_minus_it) * delta_dot / delta_f;
	  // equivalent to equation (2) of the article
	  *T_1d_it = *K_minus_it * delta;
	}

	// applying the actual stiffness
	*traction_it  = tangential_opening;
	*traction_it *= this->beta2_kappa;
	*traction_it += normal_opening;
	*traction_it *= *K_minus_it;
      }
    }

    // update precendent delta
    *delta_prec_it = delta;
  }

  delete [] memory_space;
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

INSTANTIATE_MATERIAL(MaterialCohesiveLinearFatigue);


__END_AKANTU__
