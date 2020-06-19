/**
 * @file   material_cohesive_linear_uncoupled.cc
 *
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 *
 * @date creation: Mon Jul 25 2016
 * @date last modification: Wed Feb 21 2018
 *
 * @brief  Linear irreversible cohesive law of mixed mode loading with
 * random stress definition for extrinsic type
 *
 *
 * Copyright (©) 2016-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include <algorithm>
#include <numeric>

/* -------------------------------------------------------------------------- */
#include "material_cohesive_linear_uncoupled.hh"
#include "solid_mechanics_model_cohesive.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
MaterialCohesiveLinearUncoupled<spatial_dimension>::
    MaterialCohesiveLinearUncoupled(SolidMechanicsModel & model, const ID & id)
    : MaterialCohesiveLinear<spatial_dimension>(model, id),
      delta_n_max("delta_n_max", *this), delta_t_max("delta_t_max", *this),
      damage_n("damage_n", *this), damage_t("damage_t", *this) {

  AKANTU_DEBUG_IN();

  this->registerParam(
      "roughness", R, Real(1.), _pat_parsable | _pat_readable,
      "Roughness to define coupling between mode II and mode I");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
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
template <UInt spatial_dimension>
void MaterialCohesiveLinearUncoupled<spatial_dimension>::computeTraction(
    const Array<Real> &, ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  delta_n_max.resize();
  delta_t_max.resize();
  damage_n.resize();
  damage_t.resize();

  /// define iterators
  auto traction_it =
      this->tractions(el_type, ghost_type).begin(spatial_dimension);

  auto traction_end =
      this->tractions(el_type, ghost_type).end(spatial_dimension);

  auto opening_it = this->opening(el_type, ghost_type).begin(spatial_dimension);
  auto contact_traction_it =
      this->contact_tractions(el_type, ghost_type).begin(spatial_dimension);
  auto contact_opening_it =
      this->contact_opening(el_type, ghost_type).begin(spatial_dimension);

  auto normal_it = this->normal.begin(spatial_dimension);
  auto sigma_c_it = this->sigma_c_eff(el_type, ghost_type).begin();
  auto delta_n_max_it = delta_n_max(el_type, ghost_type).begin();
  auto delta_t_max_it = delta_t_max(el_type, ghost_type).begin();
  auto delta_c_it = this->delta_c_eff(el_type, ghost_type).begin();
  auto damage_n_it = damage_n(el_type, ghost_type).begin();
  auto damage_t_it = damage_t(el_type, ghost_type).begin();

  auto insertion_stress_it =
      this->insertion_stress(el_type, ghost_type).begin(spatial_dimension);

  Vector<Real> normal_opening(spatial_dimension);
  Vector<Real> tangential_opening(spatial_dimension);

  /// loop on each quadrature point
  for (; traction_it != traction_end;
       ++traction_it, ++opening_it, ++contact_traction_it, ++contact_opening_it,
       ++normal_it, ++sigma_c_it, ++delta_n_max_it, ++delta_t_max_it,
       ++delta_c_it, ++damage_n_it, ++damage_t_it, ++insertion_stress_it) {

    Real normal_opening_norm, tangential_opening_norm;
    bool penetration;

    Real delta_c2_R2 = *delta_c_it * (*delta_c_it) / R / R;

    /// compute normal and tangential opening vectors
    normal_opening_norm = opening_it->dot(*normal_it);
    Vector<Real> normal_opening = *normal_it;
    normal_opening *= normal_opening_norm;

    //    std::cout<< "normal_opening_norm = " << normal_opening_norm
    //    <<std::endl;

    Vector<Real> tangential_opening = *opening_it;
    tangential_opening -= normal_opening;
    tangential_opening_norm = tangential_opening.norm();

    /// compute effective opening displacement
    Real delta_n =
        tangential_opening_norm * tangential_opening_norm * this->beta2_kappa2;
    Real delta_t =
        tangential_opening_norm * tangential_opening_norm * this->beta2_kappa2;

    penetration = normal_opening_norm < 0.0;
    if (this->contact_after_breaking == false &&
        Math::are_float_equal(*damage_n_it, 1.))
      penetration = false;

    if (penetration) {
      /// use penalty coefficient in case of penetration
      *contact_traction_it = normal_opening;
      *contact_traction_it *= this->penalty;
      *contact_opening_it = normal_opening;

      /// don't consider penetration contribution for delta
      //*opening_it = tangential_opening;
      normal_opening.clear();

    } else {
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
    } else {
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
    } else {
      shear_traction = tangential_opening;

      AKANTU_DEBUG_ASSERT(*delta_t_max_it != 0.,
                          "Division by zero, tolerance might be too low");

      shear_traction *= this->beta2_kappa;
      shear_traction *= *sigma_c_it / (*delta_t_max_it) * (1. - *damage_t_it);
    }

    *traction_it = normal_traction;
    *traction_it += shear_traction;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialCohesiveLinearUncoupled<spatial_dimension>::computeTangentTraction(
    const ElementType & el_type, Array<Real> & tangent_matrix,
    const Array<Real> &, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  /// define iterators
  auto tangent_it = tangent_matrix.begin(spatial_dimension, spatial_dimension);
  auto tangent_end = tangent_matrix.end(spatial_dimension, spatial_dimension);
  auto normal_it = this->normal.begin(spatial_dimension);
  auto opening_it = this->opening(el_type, ghost_type).begin(spatial_dimension);

  /// NB: delta_max_it points on delta_max_previous, i.e. the
  /// delta_max related to the solution of the previous incremental
  /// step
  auto delta_n_max_it = delta_n_max.previous(el_type, ghost_type).begin();
  auto delta_t_max_it = delta_t_max.previous(el_type, ghost_type).begin();
  auto sigma_c_it = this->sigma_c_eff(el_type, ghost_type).begin();
  auto delta_c_it = this->delta_c_eff(el_type, ghost_type).begin();
  auto damage_n_it = damage_n(el_type, ghost_type).begin();

  auto contact_opening_it =
      this->contact_opening(el_type, ghost_type).begin(spatial_dimension);

  Vector<Real> normal_opening(spatial_dimension);
  Vector<Real> tangential_opening(spatial_dimension);

  for (; tangent_it != tangent_end; ++tangent_it, ++normal_it, ++opening_it,
                                    ++sigma_c_it, ++delta_c_it,
                                    ++delta_n_max_it, ++delta_t_max_it,
                                    ++damage_n_it, ++contact_opening_it) {

    Real normal_opening_norm, tangential_opening_norm;
    bool penetration;
    Real delta_c2_R2 = *delta_c_it * (*delta_c_it) / R / R;

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

    Real delta_n =
        tangential_opening_norm * tangential_opening_norm * this->beta2_kappa2;
    Real delta_t =
        tangential_opening_norm * tangential_opening_norm * this->beta2_kappa2;

    penetration = normal_opening_norm < 0.0;
    if (this->contact_after_breaking == false &&
        Math::are_float_equal(*damage_n_it, 1.))
      penetration = false;

    Real derivative = 0; // derivative = d(t/delta)/ddelta
    Real T = 0;

    /// TANGENT STIFFNESS FOR NORMAL TRACTIONS
    Matrix<Real> n_outer_n(spatial_dimension, spatial_dimension);
    n_outer_n.outerProduct(*normal_it, *normal_it);

    if (penetration) {
      /// stiffness in compression given by the penalty parameter
      *tangent_it = n_outer_n;
      *tangent_it *= this->penalty;

      //*opening_it = tangential_opening;
      normal_opening.clear();
    } else {
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
      if (delta_n >= *delta_n_max_it) {
        if (delta_n <= *delta_c_it) {
          derivative = -(*sigma_c_it) / (delta_n * delta_n);
          T = *sigma_c_it * (1 - delta_n / *delta_c_it);
        } else {
          derivative = 0.;
          T = 0.;
        }
        // unloading-reloading
      } else if (delta_n < *delta_n_max_it) {
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

      const Vector<Real> & Delta_hat(normal_opening);
      Matrix<Real> prov(spatial_dimension, spatial_dimension);
      prov.outerProduct(Delta_hat, Delta_tilde);
      prov *= derivative / delta_n;
      prov += nn;

      Matrix<Real> prov_t = prov.transpose();

      *tangent_it = prov_t;
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
    if (delta_t >= *delta_t_max_it) {
      if (delta_t <= *delta_c_it) {
        derivative = -(*sigma_c_it) / (delta_t * delta_t);
        T = *sigma_c_it * (1 - delta_t / *delta_c_it);
      } else {
        derivative = 0.;
        T = 0.;
      }
      // unloading-reloading
    } else if (delta_t < *delta_t_max_it) {
      Real T_max = *sigma_c_it * (1 - *delta_t_max_it / *delta_c_it);
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

    Matrix<Real> prov_t = prov.transpose();

    *tangent_it += prov_t;
  }

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */

INSTANTIATE_MATERIAL(cohesive_linear_uncoupled,
                     MaterialCohesiveLinearUncoupled);

} // namespace akantu
