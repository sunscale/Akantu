/**
 * @file   material_damage_iterative.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 * @date   Sun Mar  8 12:54:30 2015
 *
 * @brief  Specialization of the class material damage to damage only one gauss
 * point at a time and propagate damage in a linear way. Max principal stress
 * criterion is used as a failure criterion.
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "material_orthotropic_damage_iterative.hh"
#include "communicator.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
MaterialOrthotropicDamageIterative<spatial_dimension>::
    MaterialOrthotropicDamageIterative(SolidMechanicsModel & model,
                                       const ID & id)
    : MaterialOrthotropicDamage<spatial_dimension>(model, id), Sc("Sc", *this),
      equivalent_stress("equivalent_stress", *this),
      stress_dir("equiv_stress_dir", *this), norm_max_equivalent_stress(0) {
  AKANTU_DEBUG_IN();

  this->registerParam("Sc", Sc, _pat_parsable, "critical stress threshold");
  this->registerParam("prescribed_dam", prescribed_dam, 0.1,
                      _pat_parsable | _pat_modifiable,
                      "increase of damage in every step");
  this->registerParam(
      "dam_threshold", dam_threshold, 0.8, _pat_parsable | _pat_modifiable,
      "damage threshold at which damage damage will be set to 1");

  this->use_previous_stress = true;
  this->use_previous_gradu = true;
  this->Sc.initialize(1);
  this->equivalent_stress.initialize(1);
  this->stress_dir.initialize(spatial_dimension * spatial_dimension);

  /// the Gauss point with the highest stress can only be of type _not_ghost
  q_max.ghost_type = _not_ghost;

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialOrthotropicDamageIterative<spatial_dimension>::
    computeNormalizedEquivalentStress(const Array<Real> & grad_u,
                                      ElementType el_type,
                                      GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  /// Vector to store eigenvalues of current stress tensor
  Vector<Real> eigenvalues(spatial_dimension);

  auto Sc_it = Sc(el_type).begin();
  auto equivalent_stress_it = equivalent_stress(el_type).begin();

  Array<Real>::const_matrix_iterator grad_u_it =
      grad_u.begin(spatial_dimension, spatial_dimension);
  Array<Real>::const_matrix_iterator grad_u_end =
      grad_u.end(spatial_dimension, spatial_dimension);

  Array<Real>::matrix_iterator stress_dir_it =
      this->stress_dir(el_type).begin(spatial_dimension, spatial_dimension);

  /// initialize matrix to store the stress for the criterion (only different in
  /// non-local)
  Matrix<Real> sigma(spatial_dimension, spatial_dimension);

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

  for (; grad_u_it != grad_u_end;
       ++Sc_it, ++equivalent_stress_it, ++stress_dir_it, ++grad_u_it) {
    sigma.zero();
    MaterialOrthotropicDamage<spatial_dimension>::computeStressOnQuad(
        *grad_u_it, sigma, 0.);

    /// rotate the tensors from the damage principal coordinate system to the CS
    /// of the computation
    if (!(Math::are_float_equal((*damage_iterator).trace(), 0))) {
      /// compute (1-D) and (1-D)^1/2
      this->computeOneMinusD(one_minus_D, *damage_iterator);
      this->computeSqrtOneMinusD(one_minus_D, sqrt_one_minus_D);

      this->rotateIntoComputationFrame(one_minus_D, one_minus_D_rotated,
                                       *damage_dir_it, rotation_tmp);

      this->rotateIntoComputationFrame(sqrt_one_minus_D,
                                       sqrt_one_minus_D_rotated, *damage_dir_it,
                                       rotation_tmp);
    } else {
      this->computeOneMinusD(one_minus_D_rotated, *damage_iterator);
      this->computeSqrtOneMinusD(one_minus_D_rotated, sqrt_one_minus_D_rotated);
    }

    computeDamageAndStressOnQuad(sigma, one_minus_D_rotated,
                                 sqrt_one_minus_D_rotated, *damage_iterator,
                                 first_term, third_term);

    /// compute the maximum principal stresses and their directions
    sigma.eig(eigenvalues, *stress_dir_it);
    *equivalent_stress_it = eigenvalues(0) / *(Sc_it);
    ++damage_dir_it;
    ++damage_iterator;
  }

  // for(;sigma_it != sigma_end; ++sigma_it,
  // 	++Sc_it, ++equivalent_stress_it, ++stress_dir_it) {
  //   /// compute the maximum principal stresses and their directions
  //   (*sigma_it).eig(eigenvalues, *stress_dir_it);
  //   *equivalent_stress_it = eigenvalues(0) / *(Sc_it);
  // }

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialOrthotropicDamageIterative<spatial_dimension>::computeAllStresses(
    GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  /// reset normalized maximum equivalent stress
  if (ghost_type == _not_ghost)
    norm_max_equivalent_stress = 0;

  MaterialOrthotropicDamage<spatial_dimension>::computeAllStresses(ghost_type);

  /// find global Gauss point with highest stress
  this->model.getMesh().getCommunicator().allReduce(
      norm_max_equivalent_stress, SynchronizerOperation::_max);

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialOrthotropicDamageIterative<spatial_dimension>::
    findMaxNormalizedEquivalentStress(ElementType el_type,
                                      GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  if (ghost_type == _not_ghost) {
    /// initialize the iterators for the equivalent stress and the damage
    const Array<Real> & e_stress = equivalent_stress(el_type);
    auto equivalent_stress_it = e_stress.begin();
    auto equivalent_stress_end = e_stress.end();
    Array<Real> & dam = this->damage(el_type);
    auto dam_it = dam.begin(this->spatial_dimension, this->spatial_dimension);
    auto damage_directions_it =
        this->damage_dir_vecs(el_type, _not_ghost)
            .begin(this->spatial_dimension, this->spatial_dimension);
    auto stress_dir_it = this->stress_dir(el_type, _not_ghost)
                             .begin(spatial_dimension, spatial_dimension);

    /// initialize the matrices for damage rotation results
    Matrix<Real> tmp(spatial_dimension, spatial_dimension);
    Matrix<Real> dam_in_computation_frame(spatial_dimension, spatial_dimension);
    Matrix<Real> dam_in_stress_frame(spatial_dimension, spatial_dimension);

    for (; equivalent_stress_it != equivalent_stress_end;
         ++equivalent_stress_it, ++dam_it, ++damage_directions_it,
         ++stress_dir_it) {
      /// check if max equivalent stress for this element type is greater than
      /// the current norm_max_eq_stress
      if (*equivalent_stress_it > norm_max_equivalent_stress &&
          (spatial_dimension * this->max_damage - (*dam_it).trace() >
           Math::getTolerance())) {

        if (Math::are_float_equal((*dam_it).trace(), 0)) {
          /// gauss point has not been damaged yet
          norm_max_equivalent_stress = *equivalent_stress_it;
          q_max.type = el_type;
          q_max.global_num = equivalent_stress_it - e_stress.begin();
        }

        else {
          /// find the damage increment on this Gauss point
          /// rotate damage into stress frame
          this->rotateIntoComputationFrame(*dam_it, dam_in_computation_frame,
                                           *damage_directions_it, tmp);
          this->rotateIntoNewFrame(dam_in_computation_frame,
                                   dam_in_stress_frame, *stress_dir_it, tmp);

          /// add damage increment
          dam_in_stress_frame(0, 0) += prescribed_dam;
          /// find new principal directions of damage
          Vector<Real> dam_eigenvalues(spatial_dimension);
          dam_in_stress_frame.eig(dam_eigenvalues);
          bool limit_reached = false;
          for (UInt i = 0; i < spatial_dimension; ++i) {
            if (dam_eigenvalues(i) + Math::getTolerance() > this->max_damage)
              limit_reached = true;
          }
          if (!limit_reached) {
            norm_max_equivalent_stress = *equivalent_stress_it;
            q_max.type = el_type;
            q_max.global_num = equivalent_stress_it - e_stress.begin();
          }
        }
      } /// end if equiv_stress > max_equiv_stress
    }   /// end loop over all gauss points of this element type
  }     // end if(_not_ghost)
  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialOrthotropicDamageIterative<spatial_dimension>::computeStress(
    ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  MaterialOrthotropicDamage<spatial_dimension>::computeStress(el_type,
                                                              ghost_type);

  auto damage_iterator =
      this->damage(el_type, ghost_type)
          .begin(this->spatial_dimension, this->spatial_dimension);
  auto damage_dir_it =
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
    this->computeOneMinusD(one_minus_D_rotated, *damage_iterator);
    this->computeSqrtOneMinusD(one_minus_D_rotated, sqrt_one_minus_D_rotated);
  }

  computeDamageAndStressOnQuad(sigma, one_minus_D_rotated,
                               sqrt_one_minus_D_rotated, *damage_iterator,
                               first_term, third_term);

  ++damage_dir_it;
  ++damage_iterator;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  computeNormalizedEquivalentStress(this->gradu(el_type, ghost_type), el_type,
                                    ghost_type);
  norm_max_equivalent_stress = 0;
  findMaxNormalizedEquivalentStress(el_type, ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
UInt MaterialOrthotropicDamageIterative<spatial_dimension>::updateDamage() {
  UInt nb_damaged_elements = 0;
  AKANTU_DEBUG_ASSERT(prescribed_dam > 0.,
                      "Your prescribed damage must be greater than zero");

  if (norm_max_equivalent_stress >= 1.) {
    AKANTU_DEBUG_IN();
    /// get the arrays and iterators for the element_type of the highest
    /// quadrature point
    ElementType el_type = q_max.type;
    UInt q_global_num = q_max.global_num;
    Array<Real> & dam = this->damage(el_type, _not_ghost);
    auto dam_it = dam.begin(this->spatial_dimension, this->spatial_dimension);
    auto damage_directions_it =
        this->damage_dir_vecs(el_type, _not_ghost)
            .begin(this->spatial_dimension, this->spatial_dimension);
    auto stress_dir_it = this->stress_dir(el_type, _not_ghost)
                             .begin(spatial_dimension, spatial_dimension);

    /// initialize the matrices for damage rotation results
    Matrix<Real> tmp(spatial_dimension, spatial_dimension);
    Matrix<Real> dam_in_computation_frame(spatial_dimension, spatial_dimension);
    Matrix<Real> dam_in_stress_frame(spatial_dimension, spatial_dimension);

    /// references to damage and directions of highest Gauss point
    Matrix<Real> q_dam = dam_it[q_global_num];
    Matrix<Real> q_dam_dir = damage_directions_it[q_global_num];
    Matrix<Real> q_stress_dir = stress_dir_it[q_global_num];

    /// increment damage
    /// find the damage increment on this Gauss point
    /// rotate damage into stress frame
    this->rotateIntoComputationFrame(q_dam, dam_in_computation_frame, q_dam_dir,
                                     tmp);
    this->rotateIntoNewFrame(dam_in_computation_frame, dam_in_stress_frame,
                             q_stress_dir, tmp);

    /// add damage increment
    dam_in_stress_frame(0, 0) += prescribed_dam;
    /// find new principal directions of damage
    Vector<Real> dam_eigenvalues(spatial_dimension);
    dam_in_stress_frame.eig(dam_eigenvalues, q_dam_dir);
    for (UInt i = 0; i < spatial_dimension; ++i) {
      q_dam(i, i) = dam_eigenvalues(i);
      if (q_dam(i, i) + Math::getTolerance() >= dam_threshold)
        q_dam(i, i) = this->max_damage;
    }
    nb_damaged_elements += 1;
  }

  this->model.getMesh().getCommunicator().allReduce(
      nb_damaged_elements, SynchronizerOperation::_sum);
  AKANTU_DEBUG_OUT();
  return nb_damaged_elements;
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialOrthotropicDamageIterative<
    spatial_dimension>::updateEnergiesAfterDamage(ElementType el_type) {
  MaterialOrthotropicDamage<spatial_dimension>::updateEnergies(el_type);
}

/* -------------------------------------------------------------------------- */
INSTANTIATE_MATERIAL(orthotropic_damage_iterative,
                     MaterialOrthotropicDamageIterative);

} // namespace akantu
