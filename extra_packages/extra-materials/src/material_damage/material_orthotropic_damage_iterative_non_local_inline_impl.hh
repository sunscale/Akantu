/**
 * @file   material_orthotropic_damage_iterative_non_local_inline_impl.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief MaterialOrthotropicDamageIterativeNonLocal inline function
 * implementation
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
} // namespace akantu

#if defined(AKANTU_DEBUG_TOOLS)
#include "aka_debug_tools.hh"
#include <string>
#endif

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
MaterialOrthotropicDamageIterativeNonLocal<spatial_dimension>::
    MaterialOrthotropicDamageIterativeNonLocal(SolidMechanicsModel & model,
                                               const ID & id)
    : Material(model, id),
      MaterialOrthotropicDamageIterativeNonLocalParent(model, id),
      grad_u_nl("grad_u non local", *this) {
  AKANTU_DEBUG_IN();
  this->is_non_local = true;
  this->grad_u_nl.initialize(spatial_dimension * spatial_dimension);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialOrthotropicDamageIterativeNonLocal<
    spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();
  this->model.getNonLocalManager().registerNonLocalVariable(
      this->gradu.getName(), grad_u_nl.getName(),
      spatial_dimension * spatial_dimension);
  MaterialOrthotropicDamageIterativeNonLocalParent::initMaterial();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialOrthotropicDamageIterativeNonLocal<
    spatial_dimension>::computeStress(ElementType type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialOrthotropicDamageIterativeNonLocal<
    spatial_dimension>::computeNonLocalStress(ElementType el_type,
                                              GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  MaterialOrthotropicDamage<spatial_dimension>::computeStress(el_type,
                                                              ghost_type);

  Array<Real>::matrix_iterator damage_iterator =
      this->damage(el_type, ghost_type)
          .begin(this->spatial_dimension, this->spatial_dimension);
  Array<Real>::matrix_iterator damage_dir_it =
      this->damage_dir_vecs(el_type, ghost_type)
          .begin(this->spatial_dimension, this->spatial_dimension);

  /// for the computation of the Cauchy stress the matrices (1-D) and
  /// (1-D)^(1/2) are needed. For the formulation see Engineering
  /// Damage Mechanics by Lemaitre and Desmorat.

  Matrix<Real> one_minus_D(this->spatial_dimension, this->spatial_dimension);
  Matrix<Real> sqrt_one_minus_D(this->spatial_dimension,
                                this->spatial_dimension);
  Matrix<Real> one_minus_D_rotated(this->spatial_dimension,
                                   this->spatial_dimension);
  Matrix<Real> sqrt_one_minus_D_rotated(this->spatial_dimension,
                                        this->spatial_dimension);
  Matrix<Real> rotation_tmp(this->spatial_dimension, this->spatial_dimension);

  /// create matrix to store the first term of the computation of the
  /// Cauchy stress
  Matrix<Real> first_term(this->spatial_dimension, this->spatial_dimension);
  Matrix<Real> third_term(this->spatial_dimension, this->spatial_dimension);

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  /// rotate the tensors from the damage principal coordinate system to the CS
  /// of the computation
  if (!(Math::are_float_equal((*damage_iterator).trace(), 0))) {
    /// compute (1-D) and (1-D)^1/2
    this->computeOneMinusD(one_minus_D, *damage_iterator);
    this->computeSqrtOneMinusD(one_minus_D, sqrt_one_minus_D);

    this->rotateIntoComputationFrame(one_minus_D, one_minus_D_rotated,
                                     *damage_dir_it, rotation_tmp);

    this->rotateIntoComputationFrame(sqrt_one_minus_D, sqrt_one_minus_D_rotated,
                                     *damage_dir_it, rotation_tmp);
  } else {
    MaterialOrthotropicDamage<spatial_dimension>::computeOneMinusD(
        one_minus_D_rotated, *damage_iterator);
    MaterialOrthotropicDamage<spatial_dimension>::computeSqrtOneMinusD(
        one_minus_D_rotated, sqrt_one_minus_D_rotated);
  }

  this->computeDamageAndStressOnQuad(sigma, one_minus_D_rotated,
                                     sqrt_one_minus_D_rotated, *damage_iterator,
                                     first_term, third_term);

  ++damage_dir_it;
  ++damage_iterator;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  this->computeNormalizedEquivalentStress(this->grad_u_nl(el_type, ghost_type),
                                          el_type, ghost_type);
  this->norm_max_equivalent_stress = 0;
  this->findMaxNormalizedEquivalentStress(el_type, ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialOrthotropicDamageIterativeNonLocal<
    spatial_dimension>::nonLocalVariableToNeighborhood() {
  this->model.getNonLocalManager().nonLocalVariableToNeighborhood(
      grad_u_nl.getName(), this->name);
}
