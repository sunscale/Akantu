/**
 * @file   custom_non_local_test_material.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sun Mar 01 2015
 * @date last modification: Mon Sep 11 2017
 *
 * @brief  Custom material to test the non local implementation
 *
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

#include "custom_non_local_test_material.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt dim>
CustomNonLocalTestMaterial<dim>::CustomNonLocalTestMaterial(
    SolidMechanicsModel & model, const ID & id)
    : MyNonLocalParent(model, id), local_damage("local_damage", *this),
      damage("damage", *this) {
  // Initialize the internal field by specifying the number of components
  this->local_damage.initialize(1);
  this->damage.initialize(1);
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void CustomNonLocalTestMaterial<dim>::registerNonLocalVariables() {
  /// register the non-local variable in the manager
  this->model.getNonLocalManager().registerNonLocalVariable(
      this->local_damage.getName(), this->damage.getName(), 1);

  this->model.getNonLocalManager()
      .getNeighborhood(this->name)
      .registerNonLocalVariable(damage.getName());
}

/* -------------------------------------------------------------------------- */
template <UInt dim> void CustomNonLocalTestMaterial<dim>::initMaterial() {
  MyNonLocalParent::initMaterial();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void CustomNonLocalTestMaterial<dim>::computeStress(ElementType el_type,
                                                    GhostType ghost_type) {
  MyNonLocalParent::computeStress(el_type, ghost_type);
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void CustomNonLocalTestMaterial<dim>::computeNonLocalStress(
    ElementType el_type, GhostType ghost_type) {
  Array<Real>::const_scalar_iterator dam =
      this->damage(el_type, ghost_type).begin();
  Array<Real>::matrix_iterator stress =
      this->stress(el_type, ghost_type).begin(dim, dim);
  Array<Real>::matrix_iterator stress_end =
      this->stress(el_type, ghost_type).end(dim, dim);

  // compute the damage and update the stresses
  for (; stress != stress_end; ++stress, ++dam) {
    *stress *= (1. - *dam);
  }
}

/* -------------------------------------------------------------------------- */
// Instantiate the material for the 3 dimensions
INSTANTIATE_MATERIAL(custom_non_local_test_material,
                     CustomNonLocalTestMaterial);
/* -------------------------------------------------------------------------- */
} // namespace akantu
