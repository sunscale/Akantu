/**
 * @file   material_non_local_inline_impl.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Aug 31 2011
 * @date last modification: Mon Jun 23 2014
 *
 * @brief  Non-local inline implementation
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
__END_AKANTU__

/* -------------------------------------------------------------------------- */
#include "aka_types.hh"
#include "integrator.hh"
#include "dumper_paraview.hh"
#include "non_local_manager.hh"
/* -------------------------------------------------------------------------- */
#include <fstream>
/* -------------------------------------------------------------------------- */

#if defined(AKANTU_DEBUG_TOOLS)
#  include "aka_debug_tools.hh"
#endif


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt DIM>
MaterialNonLocal<DIM>::MaterialNonLocal(SolidMechanicsModel & model,
					const ID & id)  :
  Material(model, id) {
  AKANTU_DEBUG_IN();

  this->model->getNonLocalManager().registerNonLocalMaterial(*this);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialNonLocal<spatial_dimension>::~MaterialNonLocal() {
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialNonLocal<spatial_dimension>::initMaterial() {
  this->insertQuadsInNeighborhoods(_not_ghost);
}

// /* -------------------------------------------------------------------------- */
// template<UInt spatial_dimension>
// void MaterialNonLocal<spatial_dimension>::updateResidual(GhostType ghost_type) {
//   AKANTU_DEBUG_IN();

//   // Update the weights for the non local variable averaging
//   if(ghost_type == _not_ghost &&
//      this->weight_func->getUpdateRate() &&
//      (this->compute_stress_calls % this->weight_func->getUpdateRate() == 0)) {
//     ElementTypeMapArray<Real> quadrature_points_coordinates("quadrature_points_coordinates", getID());
//     Mesh & mesh = this->model->getFEEngine().getMesh();
//     mesh.initElementTypeMapArray(quadrature_points_coordinates, spatial_dimension, spatial_dimension);
//     this->fem->computeQuadraturePointsCoordinates(quadrature_points_coordinates, &element_filter);
//     computeWeights(quadrature_points_coordinates);
//   }
//   if(ghost_type == _not_ghost) ++this->compute_stress_calls;

//   computeAllStresses(ghost_type);

//   computeNonLocalStresses(ghost_type);
//   assembleResidual(ghost_type);

//   AKANTU_DEBUG_OUT();
// }


// /* -------------------------------------------------------------------------- */
// template<UInt spatial_dimension>
// inline void MaterialNonLocal<spatial_dimension>::onElementsAdded(const Array<Element> & element_list) {
//   AKANTU_DEBUG_IN();
//   AKANTU_DEBUG_ERROR("This is a case not taken into account!!!");
//   AKANTU_DEBUG_OUT();
// }

// /* -------------------------------------------------------------------------- */
// template<UInt spatial_dimension>
// inline void MaterialNonLocal<spatial_dimension>::onElementsRemoved(const Array<Element> & element_list,
// 										   const ElementTypeMapArray<UInt> & new_numbering,
// 										   __attribute__((unused)) const RemovedElementsEvent & event) {
//   AKANTU_DEBUG_IN();

//   // Change the pairs in new global numbering
//   for(UInt gt = _not_ghost; gt <= _ghost; ++gt) {
//     GhostType ghost_type2 = (GhostType) gt;

//     PairList::iterator first_pair = pair_list[ghost_type2].begin();
//     PairList::iterator last_pair  = pair_list[ghost_type2].end();

//     //   Array<Real>::vector_iterator weight_it = pair_weight[ghost_type2]->begin(2);

//     for(;first_pair != last_pair; ++first_pair) {
//       QuadraturePoint & q1 = first_pair->first;
//       QuadraturePoint gq1  = this->convertToGlobalPoint(q1);
//       q1 = gq1;

//       if(new_numbering.exists(q1.type, q1.ghost_type)) {
// 	UInt q1_new_el = new_numbering(q1.type, q1.ghost_type)(gq1.element);
// 	AKANTU_DEBUG_ASSERT(q1_new_el != UInt(-1), "A local quadrature_point as been removed instead of just being renumbered");
// 	q1.element = q1_new_el;
//       }


//       QuadraturePoint & q2 = first_pair->second;
//       QuadraturePoint gq2  = this->convertToGlobalPoint(q2);
//       q2 = gq2;

//       if(new_numbering.exists(q2.type, q2.ghost_type)) {
// 	UInt q2_new_el = new_numbering(q2.type, q2.ghost_type)(gq2.element);
// 	AKANTU_DEBUG_ASSERT(q2_new_el != UInt(-1), "A local quadrature_point as been removed instead of just being renumbered");
// 	q2.element = q2_new_el;
//       }
//     }
//   }

//   // Change the material numbering
//   Material::onElementsRemoved(element_list, new_numbering, event);

//   // Change back the pairs to the new material numbering
//   for(UInt gt = _not_ghost; gt <= _ghost; ++gt) {
//     GhostType ghost_type2 = (GhostType) gt;

//     PairList::iterator first_pair = pair_list[ghost_type2].begin();
//     PairList::iterator last_pair  = pair_list[ghost_type2].end();

//     //   Array<Real>::vector_iterator weight_it = pair_weight[ghost_type2]->begin(2);

//     for(;first_pair != last_pair; ++first_pair) {
//       first_pair->first  = this->convertToLocalPoint(first_pair->first );
//       first_pair->second = this->convertToLocalPoint(first_pair->second);
//     }
//   }

//   AKANTU_DEBUG_OUT();
// }

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialNonLocal<spatial_dimension>::insertQuadsInNeighborhoods(GhostType ghost_type) {

  NonLocalManager & manager = this->model->getNonLocalManager();
  UInt spatial_dimension = this->model->getSpatialDimension();
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
	  manager.insertQuad(q, *quad, this->name);
	  ++quad;
	}
	++elem;
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialNonLocal<spatial_dimension>::updateNonLocalInternals(ElementTypeMapReal & non_local_flattened, const ID & field_id, const UInt nb_component) {

  for (ghost_type_t::iterator g = ghost_type_t::begin(); g != ghost_type_t::end(); ++g) {
    GhostType ghost_type = *g;
    /// loop over all types in the material
    typedef ElementTypeMapArray<UInt>:: type_iterator iterator;
    iterator it = this->element_filter.firstType(spatial_dimension, ghost_type, _ek_regular);
    iterator last_type = this->element_filter.lastType(spatial_dimension, ghost_type, _ek_regular);
    for(; it != last_type; ++it) {
      ElementType el_type = *it;
      Array<Real> & internal = this->getInternal<Real>(field_id)(el_type, ghost_type);
      Array<Real>::vector_iterator internal_it = internal.begin(nb_component);
      Array<Real> & internal_flat = non_local_flattened(el_type, ghost_type);
      Array<Real>::const_vector_iterator internal_flat_it = internal_flat.begin(nb_component);
      /// loop all elements for the given type
      const Array<UInt> & filter   = this->element_filter(el_type,ghost_type);
      UInt nb_elements = filter.getSize();
      UInt nb_quads = this->getFEEngine().getNbQuadraturePoints(el_type, ghost_type);
      for (UInt e = 0; e < nb_elements; ++e) {
	UInt global_el = filter(e);
	for (UInt q = 0; q < nb_quads; ++q, ++internal_it) {
	  UInt global_quad = global_el * nb_quads + q;
	  *internal_it = internal_flat_it[global_quad];
	}
      }
    }
  }  
}
