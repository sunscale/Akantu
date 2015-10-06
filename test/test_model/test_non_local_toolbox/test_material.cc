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
template<UInt spatial_dimension>
void TestMaterial<spatial_dimension>::insertQuadsInNeighborhoods(GhostType ghost_type) {

  /// this function will add all the quadrature points to the same
  /// default neighborhood instead of using one neighborhood per
  /// material
  NonLocalManager & manager = this->model->getNonLocalManager();
  InternalField<Real> quadrature_points_coordinates("quadrature_points_coordinates_tmp_nl", *this);
  quadrature_points_coordinates.initialize(spatial_dimension);

  /// intialize quadrature point object
  QuadraturePoint q;
  q.ghost_type = ghost_type;
  q.kind = _ek_regular;

  Mesh::type_iterator it = this->element_filter.firstType(spatial_dimension, ghost_type, _ek_regular);
  Mesh::type_iterator last_type = this->element_filter.lastType(spatial_dimension, ghost_type, _ek_regular);
  for(; it != last_type; ++it) {
    q.type = *it;
    const Array<UInt> & elem_filter = this->element_filter(*it, ghost_type);
    UInt nb_element  = elem_filter.getSize();
    if(nb_element) {
      UInt nb_quad = this->fem->getNbQuadraturePoints(*it, ghost_type);
      UInt nb_tot_quad = nb_quad * nb_element;
      
      Array<Real> & quads = quadrature_points_coordinates(*it, ghost_type);
      quads.resize(nb_tot_quad);

      this->model->getFEEngine().computeQuadraturePointsCoordinates(quads, *it, ghost_type, elem_filter);
      
      Array<Real>::const_vector_iterator quad = quads.begin(spatial_dimension);
      UInt * elem = elem_filter.storage();

      for (UInt e = 0; e < nb_element; ++e) {
	q.element = *elem;
	for (UInt nq = 0; nq < nb_quad; ++nq) {
	  q.num_point = nq;
	  q.global_num = q.element * nb_quad + nq;
	  manager.insertQuad(q, *quad, "test_region");
	  ++quad;
	}
	++elem;
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
// Instantiate the material for the 3 dimensions
INSTANTIATE_MATERIAL(TestMaterial);
/* -------------------------------------------------------------------------- */

__END_AKANTU__


