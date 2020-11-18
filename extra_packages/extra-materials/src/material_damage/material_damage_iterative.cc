/**
 * @file   material_damage_iterative.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
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
#include "material_damage_iterative.hh"
#include "communicator.hh"
#include "data_accessor.hh"
#include "solid_mechanics_model_RVE.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
MaterialDamageIterative<spatial_dimension>::MaterialDamageIterative(
    SolidMechanicsModel & model, const ID & id)
    : MaterialDamage<spatial_dimension>(model, id), Sc("Sc", *this),
      reduction_step("damage_step", *this),
      equivalent_stress("equivalent_stress", *this), max_reductions(0),
      norm_max_equivalent_stress(0) {
  AKANTU_DEBUG_IN();

  this->registerParam("Sc", Sc, _pat_parsable, "critical stress threshold");
  this->registerParam("prescribed_dam", prescribed_dam, 0.1,
                      _pat_parsable | _pat_modifiable, "prescribed damage");
  this->registerParam(
      "dam_threshold", dam_threshold, 0.8, _pat_parsable | _pat_modifiable,
      "damage threshold at which damage damage will be set to 1");
  this->registerParam(
      "dam_tolerance", dam_tolerance, 0.01, _pat_parsable | _pat_modifiable,
      "damage tolerance to decide if quadrature point will be damageed");
  this->registerParam("max_damage", max_damage, 0.99999,
                      _pat_parsable | _pat_modifiable, "maximum damage value");
  this->registerParam("max_reductions", max_reductions, UInt(10),
                      _pat_parsable | _pat_modifiable, "max reductions");

  this->use_previous_stress = true;
  this->use_previous_gradu = true;
  this->Sc.initialize(1);
  this->equivalent_stress.initialize(1);
  this->reduction_step.initialize(1);

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDamageIterative<spatial_dimension>::
    computeNormalizedEquivalentStress(const Array<Real> & grad_u,
                                      ElementType el_type,
                                      GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  /// Vector to store eigenvalues of current stress tensor
  Vector<Real> eigenvalues(spatial_dimension);

  auto Sc_it = Sc(el_type, ghost_type).begin();
  auto equivalent_stress_it = equivalent_stress(el_type, ghost_type).begin();

  Array<Real>::const_matrix_iterator grad_u_it =
      grad_u.begin(spatial_dimension, spatial_dimension);
  Array<Real>::const_matrix_iterator grad_u_end =
      grad_u.end(spatial_dimension, spatial_dimension);
  Real * dam = this->damage(el_type, ghost_type).storage();
  Matrix<Real> sigma(spatial_dimension, spatial_dimension);
  for (; grad_u_it != grad_u_end; ++grad_u_it) {
    sigma.zero();
    MaterialElastic<spatial_dimension>::computeStressOnQuad(*grad_u_it, sigma,
                                                            0.);
    computeDamageAndStressOnQuad(sigma, *dam);

    /// compute eigenvalues
    sigma.eig(eigenvalues);
    /// find max eigenvalue and normalize by tensile strength
    *equivalent_stress_it =
        *(std::max_element(eigenvalues.storage(),
                           eigenvalues.storage() + spatial_dimension)) /
        *(Sc_it);
    ++Sc_it;
    ++equivalent_stress_it;
    ++dam;
  }

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDamageIterative<spatial_dimension>::computeAllStresses(
    GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  /// reset normalized maximum equivalent stress
  if (ghost_type == _not_ghost)
    norm_max_equivalent_stress = 0;

  MaterialDamage<spatial_dimension>::computeAllStresses(ghost_type);

  /// find global Gauss point with highest stress
  auto rve_model = dynamic_cast<SolidMechanicsModelRVE *>(&this->model);
  if (rve_model == NULL) {
    /// is no RVE model
    const auto & comm = this->model.getMesh().getCommunicator();
    comm.allReduce(norm_max_equivalent_stress, SynchronizerOperation::_max);
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDamageIterative<spatial_dimension>::
    findMaxNormalizedEquivalentStress(ElementType el_type,
                                      GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  if (ghost_type == _not_ghost) {

    // const Array<Real> & e_stress = equivalent_stress(el_type);

    // if (e_stress.begin() != e_stress.end() ) {
    //   auto equivalent_stress_it_max =
    //   std::max_element(e_stress.begin(),e_stress.end());
    //   /// check if max equivalent stress for this element type is greater
    //   than the current norm_max_eq_stress
    //   if (*equivalent_stress_it_max > norm_max_equivalent_stress)
    // 	norm_max_equivalent_stress = *equivalent_stress_it_max;
    // }
    const Array<Real> & e_stress = equivalent_stress(el_type);
    auto equivalent_stress_it = e_stress.begin();
    auto equivalent_stress_end = e_stress.end();
    Array<Real> & dam = this->damage(el_type);
    auto dam_it = dam.begin();

    for (; equivalent_stress_it != equivalent_stress_end;
         ++equivalent_stress_it, ++dam_it) {
      /// check if max equivalent stress for this element type is greater than
      /// the current norm_max_eq_stress and if the element is not already fully
      /// damaged
      if (*equivalent_stress_it > norm_max_equivalent_stress &&
          *dam_it < max_damage) {
        norm_max_equivalent_stress = *equivalent_stress_it;
      }
    }
  }
  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDamageIterative<spatial_dimension>::computeStress(
    ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  MaterialDamage<spatial_dimension>::computeStress(el_type, ghost_type);

  Real * dam = this->damage(el_type, ghost_type).storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  computeDamageAndStressOnQuad(sigma, *dam);
  ++dam;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  computeNormalizedEquivalentStress(this->gradu(el_type, ghost_type), el_type,
                                    ghost_type);
  norm_max_equivalent_stress = 0;
  findMaxNormalizedEquivalentStress(el_type, ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
UInt MaterialDamageIterative<spatial_dimension>::updateDamage() {
  UInt nb_damaged_elements = 0;
  AKANTU_DEBUG_ASSERT(prescribed_dam > 0.,
                      "Your prescribed damage must be greater than zero");

  if (norm_max_equivalent_stress >= 1.) {

    AKANTU_DEBUG_IN();

    /// update the damage only on non-ghosts elements! Doesn't make sense to
    /// update on ghost.
    GhostType ghost_type = _not_ghost;

    for (auto && el_type : this->model.getFEEngine().getMesh().elementTypes(
             spatial_dimension, ghost_type)) {
      const Array<Real> & e_stress = equivalent_stress(el_type);
      auto equivalent_stress_it = e_stress.begin();
      auto equivalent_stress_end = e_stress.end();
      Array<Real> & dam = this->damage(el_type);
      auto dam_it = dam.begin();
      auto reduction_it = this->reduction_step(el_type, ghost_type).begin();

      for (; equivalent_stress_it != equivalent_stress_end;
           ++equivalent_stress_it, ++dam_it, ++reduction_it) {

        /// check if damage occurs
        if (*equivalent_stress_it >=
            (1 - dam_tolerance) * norm_max_equivalent_stress) {
          /// check if this element can still be damaged
          if (*reduction_it == this->max_reductions)
            continue;
          *reduction_it += 1;
          if (*reduction_it == this->max_reductions) {
            *dam_it = max_damage;
          } else {
            *dam_it += prescribed_dam;
          }
          nb_damaged_elements += 1;
        }
      }
    }
  }

  auto * rve_model = dynamic_cast<SolidMechanicsModelRVE *>(&this->model);
  if (rve_model == NULL) {
    const auto & comm = this->model.getMesh().getCommunicator();
    comm.allReduce(nb_damaged_elements, SynchronizerOperation::_sum);
  }

  AKANTU_DEBUG_OUT();
  return nb_damaged_elements;
}
/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDamageIterative<spatial_dimension>::updateEnergiesAfterDamage(
    ElementType el_type) {
  MaterialDamage<spatial_dimension>::updateEnergies(el_type);
}

/* -------------------------------------------------------------------------- */

INSTANTIATE_MATERIAL(damage_iterative, MaterialDamageIterative);

} // namespace akantu
