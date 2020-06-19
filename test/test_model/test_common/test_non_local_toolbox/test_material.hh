/**
 * @file   test_material.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Mon Sep 11 2017
 *
 * @brief  test material for the non-local neighborhood base test
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_damage.hh"
#include "material_damage_non_local.hh"

#ifndef __TEST_MATERIAL_HH__
#define __TEST_MATERIAL_HH__

using namespace akantu;

template <UInt dim>
class TestMaterial
    : public MaterialDamageNonLocal<dim, MaterialDamage<dim, MaterialElastic>> {
  /* ------------------------------------------------------------------------ */
  /* Constructor/Destructor */
  /* ------------------------------------------------------------------------ */
public:
  using Parent =
      MaterialDamageNonLocal<dim, MaterialDamage<dim, MaterialElastic>>;

  TestMaterial(SolidMechanicsModel & model, const ID & id);

  /* ------------------------------------------------------------------------ */
  /* Methods */
  /* ------------------------------------------------------------------------ */
public:
  void registerNonLocalVariables() override final;

  void computeNonLocalStress(ElementType, GhostType) override final{};

  void computeNonLocalStresses(GhostType) override final{};

protected:
  ID getNeighborhoodName() override { return "test_region"; }

  /* ------------------------------------------------------------------------ */
  /* Members */
  /* ------------------------------------------------------------------------ */
private:
  InternalField<Real> grad_u_nl;
};

#endif /* __TEST_MATERIAL_HH__ */
