/**
 * @file   test_material.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Wed Sep 23 17:16:30 2015
 *
 * @brief  Implementation of test material for the non-local neighborhood base test
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
#include "test_material.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt dim>
TestMaterial<dim>::TestMaterial(SolidMechanicsModel & model, const ID & id) :
  Material(model, id),
  MyElasticParent(model, id),
  MyNonLocalParent(model, id),
  grad_u_nl("grad_u non local", *this) {
  this->is_non_local = true;
  this->grad_u_nl.initialize(dim*dim);
 }

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void TestMaterial<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();
  this->registerNonLocalVariable(this->gradu, grad_u_nl, spatial_dimension*spatial_dimension);
  this->model->getNonLocalManager().registerNonLocalVariable(this->gradu.getName(), grad_u_nl.getName(), spatial_dimension*spatial_dimension);
  MyElasticParent::initMaterial();
  MyNonLocalParent::initMaterial();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
// Instantiate the material for the 3 dimensions
INSTANTIATE_MATERIAL(TestMaterial);
/* -------------------------------------------------------------------------- */

__END_AKANTU__


