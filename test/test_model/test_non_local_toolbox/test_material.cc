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
  MyNonLocalParent(model, id){ }

/* -------------------------------------------------------------------------- */
// template<UInt dim>
// void TestMaterial<dim>::insertQuads(NonLocalNeighborhoodBase & neighborhood) {

//   UInt spatial_dimension = this->model->getSpatialDimension();
//   InternalField<Real> quadrature_points_coordinates("quadrature_points_coordinates_tmp_nl", *this);
//   quadrature_points_coordinates.initialize(spatial_dimension);

//   GhostType ghost_type = _not_ghost;
//   QuadraturePoint q;
//   q.ghost_type = ghost_type;
//   q.kind = _ek_regular;

//   Mesh::type_iterator it        = this->element_filter.firstType(spatial_dimension, ghost_type);
//   Mesh::type_iterator last_type = this->element_filter.lastType (spatial_dimension, ghost_type);
//   for(; it != last_type; ++it) {
//     q.type = *it;
//     Array<UInt> & elem_filter = this->element_filter(*it, ghost_type);
//     UInt nb_element = elem_filter.getSize();
//     if(nb_element) {
//       UInt nb_quad = this->fem->getNbQuadraturePoints(*it, ghost_type);
//       UInt nb_tot_quad = nb_quad * nb_element;
      
//       Array<Real> & quads = quadrature_points_coordinates(*it, ghost_type);
//       quads.resize(nb_tot_quad);

//       this->model->getFEEngine().computeQuadraturePointsCoordinates(quads, *it, ghost_type, elem_filter);
      
//       Array<Real>::const_vector_iterator quad = quads.begin(spatial_dimension);
//       UInt * elem = elem_filter.storage();

//       for (UInt e = 0; e < nb_element; ++e) {
// 	q.element = *elem;
// 	for (UInt nq = 0; nq < nb_quad; ++nq) {
// 	  q.num_point = nq;
// 	  q.global_num = q.element * nb_quad + nq;
// 	  q.copyPosition(*quad);
// 	  neighborhood.insertQuad(q, *quad);
// 	  ++quad;
// 	}
// 	++elem;
//       }
//     }
//   }
// }

/* -------------------------------------------------------------------------- */
// Instantiate the material for the 3 dimensions
INSTANTIATE_MATERIAL(TestMaterial);
/* -------------------------------------------------------------------------- */

__END_AKANTU__


