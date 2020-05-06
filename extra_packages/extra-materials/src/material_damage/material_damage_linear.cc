/**
 * @file   material_damage_linear.cc
 *
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 *
 *
 * @brief  Specialization of the material class for the damage material
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "material_damage_linear.hh"
#include "solid_mechanics_model.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
MaterialDamageLinear<spatial_dimension>::MaterialDamageLinear(
    SolidMechanicsModel & model, const ID & id)
    : MaterialDamage<spatial_dimension>(model, id), K("K", *this) {
  AKANTU_DEBUG_IN();

  this->registerParam("Sigc", Sigc, 1e5, _pat_parsable, "Sigma Critique");
  this->registerParam("Gc", Gc, 2., _pat_parsable, "Gc");

  this->K.initialize(1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDamageLinear<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();
  MaterialDamage<spatial_dimension>::initMaterial();

  Epsmin = Sigc / this->E;
  Epsmax = 2 * Gc / Sigc + Epsmin;

  this->K.setDefaultValue(Epsmin);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDamageLinear<spatial_dimension>::computeStress(
    ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real * dam = this->damage(el_type, ghost_type).storage();
  Real * K = this->K(el_type, ghost_type).storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  this->computeStressOnQuad(grad_u, sigma, *dam, *K);
  ++dam;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

INSTANTIATE_MATERIAL(damage_linear, MaterialDamageLinear);

} // namespace akantu
