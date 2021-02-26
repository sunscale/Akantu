/**
 * @file   material_damage_iterative_non_local_inline_impl.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  MaterialDamageIterativeNonLocal inline function implementation
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
MaterialDamageIterativeNonLocal<spatial_dimension>::
    MaterialDamageIterativeNonLocal(SolidMechanicsModel & model, const ID & id)
    : MaterialDamageIterativeNonLocalParent(model, id),
      grad_u_nl("grad_u non local", *this) {
  AKANTU_DEBUG_IN();
  this->is_non_local = true;
  this->grad_u_nl.initialize(spatial_dimension * spatial_dimension);
  this->model.getNonLocalManager().registerNonLocalVariable(
      this->gradu.getName(), grad_u_nl.getName(),
      spatial_dimension * spatial_dimension);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDamageIterativeNonLocal<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();
  MaterialDamageIterativeNonLocalParent::initMaterial();

  this->model.getNonLocalManager().nonLocalVariableToNeighborhood(
      grad_u_nl.getName(), this->name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDamageIterativeNonLocal<spatial_dimension>::computeStress(
    ElementType /*type*/, GhostType /*ghost_type*/) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDamageIterativeNonLocal<spatial_dimension>::computeNonLocalStress(
    ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  /// compute the stress (based on the elastic law)
  MaterialDamage<spatial_dimension>::computeStress(el_type, ghost_type);

  /// multiply the stress by (1-d) to get the effective stress
  Real * dam = this->damage(el_type, ghost_type).storage();
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  this->computeDamageAndStressOnQuad(sigma, *dam);
  ++dam;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  /// compute the normalized equivalent stress
  this->computeNormalizedEquivalentStress(this->grad_u_nl(el_type, ghost_type),
                                          el_type, ghost_type);
  /// find the maximum
  this->norm_max_equivalent_stress = 0;
  this->findMaxNormalizedEquivalentStress(el_type, ghost_type);

  AKANTU_DEBUG_OUT();
}
