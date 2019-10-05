/**
 * @file   material_cohesive_exponential.cc
 *
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 * @author Seyedeh Mohadeseh Taheri Mousavi <mohadeseh.taherimousavi@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Mon Jul 09 2012
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Exponential irreversible cohesive law of mixed mode loading
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_cohesive_exponential.hh"
#include "dof_synchronizer.hh"
#include "solid_mechanics_model.hh"
#include "sparse_matrix.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
MaterialCohesiveExponential<spatial_dimension>::MaterialCohesiveExponential(
    SolidMechanicsModel & model, const ID & id)
    : MaterialCohesive(model, id) {
  AKANTU_DEBUG_IN();

  this->registerParam("beta", beta, Real(0.), _pat_parsable, "Beta parameter");

  this->registerParam("exponential_penalty", exp_penalty, true, _pat_parsable,
                      "Is contact penalty following the exponential law?");

  this->registerParam(
      "contact_tangent", contact_tangent, Real(1.0), _pat_parsable,
      "Ratio of contact tangent over the initial exponential tangent");

  // this->initInternalArray(delta_max, 1, _ek_cohesive);

  use_previous_delta_max = true;

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialCohesiveExponential<spatial_dimension>::initMaterial() {

  AKANTU_DEBUG_IN();
  MaterialCohesive::initMaterial();

  if ((exp_penalty) && (contact_tangent != 1)) {

    contact_tangent = 1;
    AKANTU_DEBUG_WARNING("The parsed paramter <contact_tangent> is forced to "
                         "1.0 when the contact penalty follows the exponential "
                         "law");
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialCohesiveExponential<spatial_dimension>::computeTraction(
    const Array<Real> & normal, ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  /// define iterators
  auto traction_it = tractions(el_type, ghost_type).begin(spatial_dimension);
  auto opening_it = opening(el_type, ghost_type).begin(spatial_dimension);
  auto normal_it = normal.begin(spatial_dimension);
  auto traction_end = tractions(el_type, ghost_type).end(spatial_dimension);
  auto delta_max_it = delta_max(el_type, ghost_type).begin();
  auto delta_max_prev_it = delta_max.previous(el_type, ghost_type).begin();

  /// compute scalars
  Real beta2 = beta * beta;

  /// loop on each quadrature point
  for (; traction_it != traction_end; ++traction_it, ++opening_it, ++normal_it,
                                      ++delta_max_it, ++delta_max_prev_it) {

    /// compute normal and tangential opening vectors
    Real normal_opening_norm = opening_it->dot(*normal_it);
    Vector<Real> normal_opening(spatial_dimension);
    normal_opening = (*normal_it);
    normal_opening *= normal_opening_norm;

    Vector<Real> tangential_opening(spatial_dimension);
    tangential_opening = *opening_it;
    tangential_opening -= normal_opening;

    Real tangential_opening_norm = tangential_opening.norm();

    /**
     * compute effective opening displacement
     * @f$ \delta = \sqrt{
     * \beta^2 \Delta_t^2 + \Delta_n^2 } @f$
     */
    Real delta = tangential_opening_norm;
    delta *= delta * beta2;
    delta += normal_opening_norm * normal_opening_norm;

    delta = sqrt(delta);

    if ((normal_opening_norm < 0) &&
        (std::abs(normal_opening_norm) > Math::getTolerance())) {

      Vector<Real> op_n(*normal_it);
      op_n *= normal_opening_norm;
      Vector<Real> delta_s(*opening_it);
      delta_s -= op_n;
      delta = tangential_opening_norm * beta;

      computeCoupledTraction(*traction_it, *normal_it, delta, delta_s,
                             *delta_max_it, *delta_max_prev_it);

      computeCompressiveTraction(*traction_it, *normal_it, normal_opening_norm,
                                 *opening_it);

    } else
      computeCoupledTraction(*traction_it, *normal_it, delta, *opening_it,
                             *delta_max_it, *delta_max_prev_it);
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialCohesiveExponential<spatial_dimension>::computeCoupledTraction(
    Vector<Real> & tract, const Vector<Real> & normal, Real delta,
    const Vector<Real> & opening, Real & delta_max_new, Real delta_max) {
  AKANTU_DEBUG_IN();

  /// full damage case
  if (std::abs(delta) < Math::getTolerance()) {
    /// set traction to zero
    tract.clear();
  } else { /// element not fully damaged
    /**
     * Compute traction loading @f$ \mathbf{T} =
     * e \sigma_c \frac{\delta}{\delta_c} e^{-\delta/ \delta_c}@f$
     */
    /**
     * Compute traction unloading @f$ \mathbf{T} =
     *  \frac{t_{max}}{\delta_{max}} \delta @f$
     */
    Real beta2 = beta * beta;
    Real normal_open_norm = opening.dot(normal);
    Vector<Real> op_n_n(spatial_dimension);
    op_n_n = normal;
    op_n_n *= (1 - beta2);
    op_n_n *= normal_open_norm;
    tract = beta2 * opening;
    tract += op_n_n;

    delta_max_new = std::max(delta_max, delta);
    tract *=
        std::exp(1.) * sigma_c * std::exp(-delta_max_new / delta_c) / delta_c;
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialCohesiveExponential<spatial_dimension>::computeCompressiveTraction(
    Vector<Real> & tract, const Vector<Real> & normal, Real delta_n,
    __attribute__((unused)) const Vector<Real> & opening) {
  Vector<Real> temp_tract(normal);

  if (exp_penalty) {
    temp_tract *= delta_n * std::exp(1) * sigma_c *
                  std::exp(-delta_n / delta_c) / delta_c;
  } else {
    Real initial_tg =
        contact_tangent * std::exp(1.) * sigma_c * delta_n / delta_c;
    temp_tract *= initial_tg;
  }

  tract += temp_tract;
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialCohesiveExponential<spatial_dimension>::computeTangentTraction(
    const ElementType & el_type, Array<Real> & tangent_matrix,
    const Array<Real> & normal, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto tangent_it = tangent_matrix.begin(spatial_dimension, spatial_dimension);
  auto tangent_end = tangent_matrix.end(spatial_dimension, spatial_dimension);
  auto normal_it = normal.begin(spatial_dimension);
  auto opening_it = opening(el_type, ghost_type).begin(spatial_dimension);
  auto delta_max_it = delta_max.previous(el_type, ghost_type).begin();

  Real beta2 = beta * beta;

  /**
   * compute tangent matrix  @f$ \frac{\partial \mathbf{t}}
   * {\partial \delta} = \hat{\mathbf{t}} \otimes
   * \frac{\partial (t/\delta)}{\partial \delta}
   * \frac{\hat{\mathbf{t}}}{\delta}+ \frac{t}{\delta}  [ \beta^2 \mathbf{I} +
   * (1-\beta^2) (\mathbf{n} \otimes \mathbf{n})] @f$
   **/

  /**
   * In which @f$
   *  \frac{\partial(t/ \delta)}{\partial \delta} =
   * \left\{\begin{array} {l l}
   *  -e  \frac{\sigma_c}{\delta_c^2  }e^{-\delta  /  \delta_c} &  \quad  if
   *  \delta \geq \delta_{max} \\
   *  0 & \quad if \delta < \delta_{max}, \delta_n > 0
   *  \end{array}\right. @f$
   **/

  for (; tangent_it != tangent_end;
       ++tangent_it, ++normal_it, ++opening_it, ++delta_max_it) {

    Real normal_opening_norm = opening_it->dot(*normal_it);
    Vector<Real> normal_opening(spatial_dimension);
    normal_opening = (*normal_it);
    normal_opening *= normal_opening_norm;

    Vector<Real> tangential_opening(spatial_dimension);
    tangential_opening = *opening_it;
    tangential_opening -= normal_opening;

    Real tangential_opening_norm = tangential_opening.norm();

    Real delta = tangential_opening_norm;
    delta *= delta * beta2;
    delta += normal_opening_norm * normal_opening_norm;
    delta = sqrt(delta);

    if ((normal_opening_norm < 0) &&
        (std::abs(normal_opening_norm) > Math::getTolerance())) {

      Vector<Real> op_n(*normal_it);
      op_n *= normal_opening_norm;
      Vector<Real> delta_s(*opening_it);
      delta_s -= op_n;
      delta = tangential_opening_norm * beta;

      computeCoupledTangent(*tangent_it, *normal_it, delta, delta_s,
                            *delta_max_it);

      computeCompressivePenalty(*tangent_it, *normal_it, normal_opening_norm);

    } else
      computeCoupledTangent(*tangent_it, *normal_it, delta, *opening_it,
                            *delta_max_it);
  }

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialCohesiveExponential<spatial_dimension>::computeCoupledTangent(
    Matrix<Real> & tangent, const Vector<Real> & normal, Real delta,
    const Vector<Real> & opening, Real) {
  AKANTU_DEBUG_IN();

  Real beta2 = beta * beta;
  Matrix<Real> J(spatial_dimension, spatial_dimension);
  J.eye(beta2);

  if (std::abs(delta) < Math::getTolerance()) {
    delta = Math::getTolerance();
  }

  Real opening_normal;
  opening_normal = opening.dot(normal);

  Vector<Real> delta_e(normal);
  delta_e *= opening_normal;
  delta_e *= (1. - beta2);
  delta_e += (beta2 * opening);

  Real exponent = std::exp(1. - delta / delta_c) * sigma_c / delta_c;

  Matrix<Real> first_term(spatial_dimension, spatial_dimension);
  first_term.outerProduct(normal, normal);
  first_term *= (1. - beta2);
  first_term += J;

  Matrix<Real> second_term(spatial_dimension, spatial_dimension);
  second_term.outerProduct(delta_e, delta_e);
  second_term /= delta;
  second_term /= delta_c;

  Matrix<Real> diff(first_term);
  diff -= second_term;

  tangent = diff;
  tangent *= exponent;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialCohesiveExponential<spatial_dimension>::computeCompressivePenalty(
    Matrix<Real> & tangent, const Vector<Real> & normal, Real delta_n) {

  if (!exp_penalty)
    delta_n = 0.;

  Matrix<Real> n_outer_n(spatial_dimension, spatial_dimension);
  n_outer_n.outerProduct(normal, normal);

  Real normal_tg = contact_tangent * std::exp(1.) * sigma_c *
                   std::exp(-delta_n / delta_c) * (1. - delta_n / delta_c) /
                   delta_c;

  n_outer_n *= normal_tg;

  tangent += n_outer_n;
}

INSTANTIATE_MATERIAL(cohesive_exponential, MaterialCohesiveExponential);

} // namespace akantu
