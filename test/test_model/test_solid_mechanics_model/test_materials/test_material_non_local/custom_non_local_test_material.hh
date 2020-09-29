/**
 * @file   custom_non_local_test_material.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Aug 23 2012
 * @date last modification: Mon Sep 11 2017
 *
 * @brief  Custom material to test the non local implementation
 *
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
#include "material_elastic.hh"
#include "material_non_local.hh"
/* -------------------------------------------------------------------------- */

#ifndef CUSTOM_NON_LOCAL_TEST_MATERIAL_HH_
#define CUSTOM_NON_LOCAL_TEST_MATERIAL_HH_

namespace akantu {

template <UInt dim>
class CustomNonLocalTestMaterial
    : public MaterialNonLocal<dim, MaterialElastic<dim>> {
public:
  using MyNonLocalParent = MaterialNonLocal<dim, MaterialElastic<dim>>;

  CustomNonLocalTestMaterial(SolidMechanicsModel & model, const ID & id);

  /* ------------------------------------------------------------------------ */
  void initMaterial() override;

  void computeNonLocalStress(ElementType el_type, GhostType ghost_type);
  void computeStress(ElementType el_type, GhostType ghost_type) override;

protected:
  void registerNonLocalVariables() override;

  /* ------------------------------------------------------------------------ */
  void computeNonLocalStresses(GhostType ghost_type) override {
    AKANTU_DEBUG_IN();

    for (auto & type : this->element_filter.elementTypes(dim, ghost_type)) {
      computeNonLocalStress(type, ghost_type);
    }

    AKANTU_DEBUG_OUT();
  }

public:
  void setDamage(Real dam) { this->local_damage.setDefaultValue(dam); }

protected:
  InternalField<Real> local_damage;
  InternalField<Real> damage;
};

} // namespace akantu

#endif /* CUSTOM_NON_LOCAL_TEST_MATERIAL_HH_ */
