/**
 * @file   material_cohesive_linear_inline_impl.cc
 *
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Apr 22 2015
 * @date last modification: Wed Feb 21 2018
 *
 * @brief  Inline functions of the MaterialCohesiveLinear
 *
 * @section LICENSE
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_cohesive_linear.hh"
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_MATERIAL_COHESIVE_LINEAR_INLINE_IMPL_CC__
#define __AKANTU_MATERIAL_COHESIVE_LINEAR_INLINE_IMPL_CC__

/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt dim>
inline Real MaterialCohesiveLinear<dim>::computeEffectiveNorm(
    const Matrix<Real> & stress, const Vector<Real> & normal,
    const Vector<Real> & tangent, Vector<Real> & normal_traction) const {
  normal_traction.mul<false>(stress, normal);

  Real normal_contrib = normal_traction.dot(normal);

  /// in 3D tangential components must be summed
  Real tangent_contrib = 0;

  if (dim == 2) {
    Real tangent_contrib_tmp = normal_traction.dot(tangent);
    tangent_contrib += tangent_contrib_tmp * tangent_contrib_tmp;
  } else if (dim == 3) {
    for (UInt s = 0; s < dim - 1; ++s) {
      const Vector<Real> tangent_v(tangent.storage() + s * dim, dim);
      Real tangent_contrib_tmp = normal_traction.dot(tangent_v);
      tangent_contrib += tangent_contrib_tmp * tangent_contrib_tmp;
    }
  }

  tangent_contrib = std::sqrt(tangent_contrib);
  normal_contrib = std::max(Real(0.), normal_contrib);

  return std::sqrt(normal_contrib * normal_contrib +
                   tangent_contrib * tangent_contrib * beta2_inv);
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
inline void MaterialCohesiveLinear<dim>::computeTractionOnQuad(
    Vector<Real> & traction, Vector<Real> & opening,
    const Vector<Real> & normal, Real & delta_max, const Real & delta_c,
    const Vector<Real> & insertion_stress, const Real & sigma_c,
    Vector<Real> & normal_opening, Vector<Real> & tangential_opening,
    Real & normal_opening_norm, Real & tangential_opening_norm, Real & damage,
    bool & penetration, Vector<Real> & contact_traction,
    Vector<Real> & contact_opening) {

  /// compute normal and tangential opening vectors
  normal_opening_norm = opening.dot(normal);
  normal_opening = normal;
  normal_opening *= normal_opening_norm;

  tangential_opening = opening;
  tangential_opening -= normal_opening;
  tangential_opening_norm = tangential_opening.norm();

  /**
   * compute effective opening displacement
   * @f$ \delta = \sqrt{
   * \frac{\beta^2}{\kappa^2} \Delta_t^2 + \Delta_n^2 } @f$
   */
  Real delta =
      tangential_opening_norm * tangential_opening_norm * this->beta2_kappa2;

  penetration = normal_opening_norm / delta_c < -Math::getTolerance();
  // penetration = normal_opening_norm < 0.;
  if (this->contact_after_breaking == false &&
      Math::are_float_equal(damage, 1.))
    penetration = false;

  if (penetration) {
    /// use penalty coefficient in case of penetration
    contact_traction = normal_opening;
    contact_traction *= this->penalty;
    contact_opening = normal_opening;

    /// don't consider penetration contribution for delta
    opening = tangential_opening;
    normal_opening.clear();
  } else {
    delta += normal_opening_norm * normal_opening_norm;
    contact_traction.clear();
    contact_opening.clear();
  }

  delta = std::sqrt(delta);

  /// update maximum displacement and damage
  delta_max = std::max(delta_max, delta);
  damage = std::min(delta_max / delta_c, Real(1.));

  /**
   * Compute traction @f$ \mathbf{T} = \left(
   * \frac{\beta^2}{\kappa} \Delta_t \mathbf{t} + \Delta_n
   * \mathbf{n} \right) \frac{\sigma_c}{\delta} \left( 1-
   * \frac{\delta}{\delta_c} \right)@f$
   */

  if (Math::are_float_equal(damage, 1.))
    traction.clear();
  else if (Math::are_float_equal(damage, 0.)) {
    if (penetration)
      traction.clear();
    else
      traction = insertion_stress;
  } else {
    traction = tangential_opening;
    traction *= this->beta2_kappa;
    traction += normal_opening;

    AKANTU_DEBUG_ASSERT(delta_max != 0.,
                        "Division by zero, tolerance might be too low");

    traction *= sigma_c / delta_max * (1. - damage);
  }
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
inline void MaterialCohesiveLinear<dim>::computeTangentTractionOnQuad(
    Matrix<Real> & tangent, Real & delta_max, const Real & delta_c,
    const Real & sigma_c, Vector<Real> & opening, const Vector<Real> & normal,
    Vector<Real> & normal_opening, Vector<Real> & tangential_opening,
    Real & normal_opening_norm, Real & tangential_opening_norm, Real & damage,
    bool & penetration, Vector<Real> & contact_opening) {

  /**
   * During the update of the residual the interpenetrations are
   * stored in the array "contact_opening", therefore, in the case
   * of penetration, in the array "opening" there are only the
   * tangential components.
   */
  opening += contact_opening;

  /// compute normal and tangential opening vectors
  normal_opening_norm = opening.dot(normal);
  normal_opening = normal;
  normal_opening *= normal_opening_norm;

  tangential_opening = opening;
  tangential_opening -= normal_opening;
  tangential_opening_norm = tangential_opening.norm();

  Real delta =
      tangential_opening_norm * tangential_opening_norm * this->beta2_kappa2;

  penetration = normal_opening_norm < 0.0;
  if (this->contact_after_breaking == false &&
      Math::are_float_equal(damage, 1.))
    penetration = false;

  Real derivative = 0; // derivative = d(t/delta)/ddelta
  Real t = 0;

  Matrix<Real> n_outer_n(spatial_dimension, spatial_dimension);
  n_outer_n.outerProduct(normal, normal);

  if (penetration) {
    /// stiffness in compression given by the penalty parameter
    tangent += n_outer_n;
    tangent *= penalty;

    opening = tangential_opening;
    normal_opening_norm = opening.dot(normal);
    normal_opening = normal;
    normal_opening *= normal_opening_norm;
  } else {
    delta += normal_opening_norm * normal_opening_norm;
  }

  delta = std::sqrt(delta);

  /**
   * Delta has to be different from 0 to have finite values of
   * tangential stiffness.  At the element insertion, delta =
   * 0. Therefore, a fictictious value is defined, for the
   * evaluation of the first value of K.
   */
  if (delta < Math::getTolerance())
    delta = delta_c / 1000.;

  if (delta >= delta_max) {
    if (delta <= delta_c) {
      derivative = -sigma_c / (delta * delta);
      t = sigma_c * (1 - delta / delta_c);
    } else {
      derivative = 0.;
      t = 0.;
    }
  } else if (delta < delta_max) {
    Real tmax = sigma_c * (1 - delta_max / delta_c);
    t = tmax / delta_max * delta;
  }

  /// computation of the derivative of the constitutive law (dT/ddelta)
  Matrix<Real> I(spatial_dimension, spatial_dimension);
  I.eye(this->beta2_kappa);

  Matrix<Real> nn(n_outer_n);
  nn *= (1. - this->beta2_kappa);
  nn += I;
  nn *= t / delta;

  Vector<Real> t_tilde(normal_opening);
  t_tilde *= (1. - this->beta2_kappa2);

  Vector<Real> mm(opening);
  mm *= this->beta2_kappa2;
  t_tilde += mm;

  Vector<Real> t_hat(normal_opening);
  t_hat += this->beta2_kappa * tangential_opening;

  Matrix<Real> prov(spatial_dimension, spatial_dimension);
  prov.outerProduct(t_hat, t_tilde);
  prov *= derivative / delta;
  prov += nn;

  Matrix<Real> prov_t = prov.transpose();
  tangent += prov_t;
}

/* -------------------------------------------------------------------------- */
} // namespace akantu

/* -------------------------------------------------------------------------- */
#endif //__AKANTU_MATERIAL_COHESIVE_LINEAR_INLINE_IMPL_CC__
