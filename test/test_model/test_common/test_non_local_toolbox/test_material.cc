/**
 * @file   test_material.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 * @date creation: Sat Sep 26 2015
 * @date last modification: Tue Sep 19 2017
 *
 * @brief  Implementation of test material for the non-local neighborhood base
 * test
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
#include "test_material.hh"

/* -------------------------------------------------------------------------- */
template <UInt dim>
TestMaterial<dim>::TestMaterial(SolidMechanicsModel & model, const ID & id)
    : Parent(model, id), grad_u_nl("grad_u non local", *this) {
  this->is_non_local = true;
  this->grad_u_nl.initialize(dim * dim);
}

/* -------------------------------------------------------------------------- */
template <UInt dim> void TestMaterial<dim>::registerNonLocalVariables() {
  this->model.getNonLocalManager().registerNonLocalVariable(
      this->gradu.getName(), grad_u_nl.getName(), dim * dim);

  this->model.getNonLocalManager()
      .getNeighborhood(this->getNeighborhoodName())
      .registerNonLocalVariable(grad_u_nl.getName());
}

/* -------------------------------------------------------------------------- */
// Instantiate the material for the 3 dimensions
INSTANTIATE_MATERIAL(test_material, TestMaterial);
/* -------------------------------------------------------------------------- */
