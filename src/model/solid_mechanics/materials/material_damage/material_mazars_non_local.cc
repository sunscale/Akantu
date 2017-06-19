/**
 * @file   material_mazars_non_local.cc
 *
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Thu Oct 15 2015
 *
 * @brief  Specialization of the material class for the non-local mazars
 * material
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "material_mazars_non_local.hh"
#include "solid_mechanics_model.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
MaterialMazarsNonLocal<spatial_dimension>::MaterialMazarsNonLocal(
    SolidMechanicsModel & model, const ID & id)
    : Material(model, id), MaterialMazars<spatial_dimension>(model, id),
      MaterialNonLocalParent(model, id), Ehat("epsilon_equ", *this) {
  AKANTU_DEBUG_IN();

  this->is_non_local = true;
  this->Ehat.initialize(1);

  this->registerParam("average_on_damage", this->damage_in_compute_stress,
                      false, _pat_parsable | _pat_modifiable,
                      "Is D the non local variable");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialMazarsNonLocal<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();
  MaterialMazars<spatial_dimension>::initMaterial();
  MaterialNonLocalParent::initMaterial();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialMazarsNonLocal<spatial_dimension>::computeStress(
    ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real * damage = this->damage(el_type, ghost_type).storage();
  Real * epsilon_equ = this->Ehat(el_type, ghost_type).storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  MaterialMazars<spatial_dimension>::computeStressOnQuad(grad_u, sigma, *damage,
                                                         *epsilon_equ);
  ++damage;
  ++epsilon_equ;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialMazarsNonLocal<spatial_dimension>::computeNonLocalStresses(
    __attribute__((unused)) GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  InternalField<Real> nl_var("Non local variable", *this);
  nl_var.initialize(1);

  // if(this->damage_in_compute_stress)
  //   this->weightedAvergageOnNeighbours(this->damage, nl_var, 1);
  // else
  //   this->weightedAvergageOnNeighbours(this->Ehat, nl_var, 1);

  // Mesh::type_iterator it =
  // this->model->getFEEngine().getMesh().firstType(spatial_dimension,
  // ghost_type);
  // Mesh::type_iterator last_type =
  // this->model->getFEEngine().getMesh().lastType(spatial_dimension,
  // ghost_type);
  // for(; it != last_type; ++it) {
  //   this->computeNonLocalStress(nl_var(*it, ghost_type), *it, ghost_type);
  // }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialMazarsNonLocal<spatial_dimension>::computeNonLocalStress(
    Array<Real> & non_loc_var, ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real * damage;
  Real * epsilon_equ;
  if (this->damage_in_compute_stress) {
    damage = non_loc_var.storage();
    epsilon_equ = this->Ehat(el_type, ghost_type).storage();
  } else {
    damage = this->damage(el_type, ghost_type).storage();
    epsilon_equ = non_loc_var.storage();
  }

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);
  this->computeDamageAndStressOnQuad(grad_u, sigma, *damage, *epsilon_equ);

  ++damage;
  ++epsilon_equ;
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialMazarsNonLocal<
    spatial_dimension>::nonLocalVariableToNeighborhood() {}

INSTANTIATE_MATERIAL(MaterialMazarsNonLocal);

} // akantu
