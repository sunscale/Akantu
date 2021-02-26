/**
 * @file   material_orthotropic_damage_tmpl.hh
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Sun Mar  8 12:54:30 2015
 *
 * @brief Specialization of the material class for the orthotropic
 * damage material
 *
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
#include "material_orthotropic_damage.hh"
#include "solid_mechanics_model.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension, template <UInt> class Parent>
MaterialOrthotropicDamage<spatial_dimension, Parent>::MaterialOrthotropicDamage(
    SolidMechanicsModel & model, const ID & id)
    : Parent<spatial_dimension>(model, id), damage("damage", *this),
      dissipated_energy("damage dissipated energy", *this),
      int_sigma("integral of sigma", *this),
      damage_dir_vecs("damage_principal_directions", *this) {
  AKANTU_DEBUG_IN();

  this->registerParam("eta", eta, 2., _pat_parsable | _pat_modifiable,
                      "Damage sensitivity parameter");
  this->registerParam("max_damage", max_damage, 0.99999,
                      _pat_parsable | _pat_modifiable, "maximum damage value");

  this->is_non_local = false;
  this->use_previous_stress = true;
  this->use_previous_gradu = true;

  /// use second order tensor for description of damage state
  this->damage.initialize(spatial_dimension * spatial_dimension);
  this->dissipated_energy.initialize(1);
  this->int_sigma.initialize(1);
  this->damage_dir_vecs.initialize(spatial_dimension * spatial_dimension);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension, template <UInt> class Parent>
void MaterialOrthotropicDamage<spatial_dimension, Parent>::initMaterial() {
  AKANTU_DEBUG_IN();
  Parent<spatial_dimension>::initMaterial();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Compute the dissipated energy in  each element by a trapezoidal approximation
 * of
 * @f$ Ed = \int_0^{\epsilon}\sigma(\omega)d\omega -
 * \frac{1}{2}\sigma:\epsilon@f$
 */
template <UInt spatial_dimension, template <UInt> class Parent>
void MaterialOrthotropicDamage<spatial_dimension, Parent>::updateEnergies(
    ElementType el_type) {
  Parent<spatial_dimension>::updateEnergies(el_type);

  this->computePotentialEnergy(el_type);

  auto epsilon_p =
      this->gradu.previous(el_type).begin(spatial_dimension, spatial_dimension);
  auto sigma_p = this->stress.previous(el_type).begin(spatial_dimension,
                                                      spatial_dimension);

  auto epot = this->potential_energy(el_type).begin();

  auto ints = this->int_sigma(el_type).begin();
  auto ed = this->dissipated_energy(el_type).begin();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, _not_ghost);

  Matrix<Real> delta_gradu_it(grad_u);
  delta_gradu_it -= *epsilon_p;

  Matrix<Real> sigma_h(sigma);
  sigma_h += *sigma_p;

  Real dint = .5 * sigma_h.doubleDot(delta_gradu_it);

  *ints += dint;
  *ed = *ints - *epot;

  ++epsilon_p;
  ++sigma_p;
  ++epot;
  ++ints;
  ++ed;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension, template <UInt> class Parent>
void MaterialOrthotropicDamage<spatial_dimension, Parent>::computeTangentModuli(
    ElementType el_type, Array<Real> & tangent_matrix,
    GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Parent<spatial_dimension>::computeTangentModuli(el_type, tangent_matrix,
                                                  ghost_type);

  /// get the damage array for current element type
  Array<Real> & dam = this->damage(el_type);
  auto dam_it = dam.begin(this->spatial_dimension, this->spatial_dimension);

  /// get the directions of damage for the current element type
  Array<Real> & dam_dirs = this->damage_dir_vecs(el_type);
  auto damage_directions_it =
      dam_dirs.begin(this->spatial_dimension, this->spatial_dimension);

  /// for the computation of the Cauchy stress the matrices (1-D) and
  /// (1-D)^(1/2) are needed. For the formulation see Engineering
  /// Damage Mechanics by Lemaitre and Desmorat.
  Matrix<Real> one_minus_D(this->spatial_dimension, this->spatial_dimension);
  Matrix<Real> sqrt_one_minus_D(this->spatial_dimension,
                                this->spatial_dimension);
  Matrix<Real> one_minus_D_rot(spatial_dimension, spatial_dimension);
  Matrix<Real> sqrt_one_minus_D_rot(spatial_dimension, spatial_dimension);

  Matrix<Real> rotation_tmp(spatial_dimension, spatial_dimension);

  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_matrix);

  if (!(Math::are_float_equal((*dam_it).trace(), 0)))
    computeTangentModuliOnQuad(tangent, tangent, *dam_it, *damage_directions_it,
                               one_minus_D, sqrt_one_minus_D, one_minus_D_rot,
                               sqrt_one_minus_D_rot, rotation_tmp);

  ++dam_it;
  ++damage_directions_it;

  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension, template <UInt> class Parent>
void MaterialOrthotropicDamage<spatial_dimension, Parent>::
    computeTangentModuliOnQuad(Matrix<Real> & tangent, const Matrix<Real> C,
                               const Matrix<Real> & dam,
                               const Matrix<Real> & dam_directions,
                               Matrix<Real> & O_prime, Matrix<Real> & S_prime,
                               Matrix<Real> & O, Matrix<Real> & S,
                               Matrix<Real> & rotation_tmp) {

  /// effect of damage on stiffness matrix: See Ragueneau et al. 2008, p. 423,
  /// ep. 7
  Real trace_D = dam.trace();
  this->computeOneMinusD(O_prime, dam);
  this->computeSqrtOneMinusD(O_prime, S_prime);
  this->rotateIntoComputationFrame(O_prime, O, dam_directions, rotation_tmp);
  this->rotateIntoComputationFrame(S_prime, S, dam_directions, rotation_tmp);

  /// compute new stiffness matrix in damage coordinate system
  if (spatial_dimension == 1)
    tangent *= (1 - dam(0, 0));

  if (spatial_dimension == 2) {
    Real min_val =
        std::min((this->eta / spatial_dimension * trace_D), this->max_damage);

    /// first row
    tangent(0, 0) =
        (C(0, 0) * S(0, 0) * S(0, 0) + C(1, 0) * S(0, 1) * S(0, 1) -
         (min_val / 2. - 1. / 2) * (C(0, 0) + C(1, 0)) +
         (O(0, 0) * (C(0, 0) * O(0, 0) + C(1, 0) * O(1, 1))) / (trace_D - 2.));

    tangent(0, 1) =
        (C(0, 1) * S(0, 0) * S(0, 0) + C(1, 1) * S(0, 1) * S(0, 1) -
         (min_val / 2. - 1. / 2) * (C(0, 1) + C(1, 1)) +
         (O(0, 0) * (C(0, 1) * O(0, 0) + C(1, 1) * O(1, 1))) / (trace_D - 2.));

    tangent(0, 2) = (2. * C(2, 2) * S(0, 0) * S(0, 1) +
                     (2. * C(2, 2) * O(0, 0) * O(0, 1)) / (trace_D - 2.));

    /// second row
    tangent(1, 0) =
        (C(0, 0) * S(0, 1) * S(0, 1) + C(1, 0) * S(1, 1) * S(1, 1) -
         (min_val / 2. - 1. / 2) * (C(0, 0) + C(1, 0)) +
         (O(1, 1) * (C(0, 0) * O(0, 0) + C(1, 0) * O(1, 1))) / (trace_D - 2.));

    tangent(1, 1) =
        (C(0, 1) * S(0, 1) * S(0, 1) + C(1, 1) * S(1, 1) * S(1, 1) -
         (min_val / 2. - 1. / 2) * (C(0, 1) + C(1, 1)) +
         (O(1, 1) * (C(0, 1) * O(0, 0) + C(1, 1) * O(1, 1))) / (trace_D - 2.));

    tangent(1, 2) = (2. * C(2, 2) * S(0, 1) * S(1, 1) +
                     (2. * C(2, 2) * O(0, 1) * O(1, 1)) / (trace_D - 2.));

    /// third row
    tangent(2, 0) =
        ((O(0, 1) * (C(0, 0) * O(0, 0) + C(1, 0) * O(1, 1))) / (trace_D - 2.) +
         C(0, 0) * S(0, 0) * S(0, 1) + C(1, 0) * S(0, 1) * S(1, 1));

    tangent(2, 1) =
        ((O(0, 1) * (C(0, 1) * O(0, 0) + C(1, 1) * O(1, 1))) / (trace_D - 2.) +
         C(0, 1) * S(0, 0) * S(0, 1) + C(1, 1) * S(0, 1) * S(1, 1));
    tangent(2, 2) = ((2. * C(2, 2) * O(0, 1) * O(0, 1)) / (trace_D - 2.) +
                     C(2, 2) * S(0, 1) * S(0, 1) + C(2, 2) * S(0, 0) * S(1, 1));
  }

  if (spatial_dimension == 3) {
  }
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension, template <UInt> class Parent>
inline void MaterialOrthotropicDamage<spatial_dimension, Parent>::
    computeDamageAndStressOnQuad(Matrix<Real> & sigma,
                                 Matrix<Real> & one_minus_D,
                                 Matrix<Real> & sqrt_one_minus_D,
                                 Matrix<Real> & damage,
                                 Matrix<Real> & first_term,
                                 Matrix<Real> & third_term) {
  /// Definition of Cauchy stress based on second order damage tensor:
  /// "Anisotropic damage modelling of biaxial behaviour and rupture
  /// of concrete strucutres", Ragueneau et al., 2008, Eq. 7
  first_term.mul<false, false>(sqrt_one_minus_D, sigma);
  first_term *= sqrt_one_minus_D;

  Real second_term = 0;
  for (UInt i = 0; i < this->spatial_dimension; ++i) {
    for (UInt j = 0; j < this->spatial_dimension; ++j)
      second_term += sigma(i, j) * one_minus_D(i, j);
  }

  second_term /= (this->spatial_dimension - damage.trace());

  one_minus_D *= second_term;

  third_term.eye(1. / this->spatial_dimension * sigma.trace() *
                 (1 - eta / (this->spatial_dimension) * damage.trace()));

  sigma.copy(first_term);
  sigma -= one_minus_D;
  sigma += third_term;
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension, template <UInt> class Parent>
inline void MaterialOrthotropicDamage<spatial_dimension, Parent>::
    rotateIntoComputationFrame(const Matrix<Real> & to_rotate,
                               Matrix<Real> & rotated,
                               const Matrix<Real> & damage_directions,
                               Matrix<Real> & rotation_tmp) {
  rotation_tmp.mul<false, true>(to_rotate, damage_directions);
  rotated.mul<false, false>(damage_directions, rotation_tmp);
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension, template <UInt> class Parent>
inline void
MaterialOrthotropicDamage<spatial_dimension, Parent>::rotateIntoNewFrame(
    const Matrix<Real> & to_rotate, Matrix<Real> & rotated,
    const Matrix<Real> & damage_directions, Matrix<Real> & rotation_tmp) {
  rotation_tmp.mul<false, false>(to_rotate, damage_directions);
  rotated.mul<true, false>(damage_directions, rotation_tmp);
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension, template <UInt> class Parent>
inline void
MaterialOrthotropicDamage<spatial_dimension, Parent>::computeOneMinusD(
    Matrix<Real> & one_minus_D, const Matrix<Real> & damage) {
  /// compute one minus
  one_minus_D.eye();
  one_minus_D -= damage;
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension, template <UInt> class Parent>
inline void
MaterialOrthotropicDamage<spatial_dimension, Parent>::computeSqrtOneMinusD(
    const Matrix<Real> & one_minus_D, Matrix<Real> & sqrt_one_minus_D) {

/// To compute (1-D)^1/2 we need to check that we are in the
/// principal coordinate system of the damage
#ifndef AKANTU_NDEBUG
  for (UInt i = 0; i < this->spatial_dimension; ++i) {
    for (UInt j = 0; j < this->spatial_dimension; ++j) {
      if (i != j)
        AKANTU_DEBUG_ASSERT(Math::are_float_equal(one_minus_D(i, j), 0),
                            "The damage tensor has off-diagonal parts");
    }
  }
#endif // AKANTU_NDEBUG

  /// compute (1-D)^1/2
  sqrt_one_minus_D.copy(one_minus_D);
  for (UInt i = 0; i < this->spatial_dimension; ++i)
    sqrt_one_minus_D(i, i) = std::sqrt(sqrt_one_minus_D(i, i));
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
