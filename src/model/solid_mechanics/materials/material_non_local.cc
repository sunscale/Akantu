/**
 * @file   material_non_local.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Thu Oct  8 15:12:27 2015
 *
 * @brief  Implementation of material non-local
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
#include "material_non_local.hh"
/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__
/* -------------------------------------------------------------------------- */
template<UInt DIM>
MaterialNonLocal<DIM>::MaterialNonLocal(SolidMechanicsModel & model,
					const ID & id)  :
  Material(model, id) {
  AKANTU_DEBUG_IN();

  NonLocalManager & manager = this->model->getNonLocalManager();
  manager.registerNonLocalMaterial(*this);
 
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialNonLocal<spatial_dimension>::~MaterialNonLocal() {
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialNonLocal<spatial_dimension>::initMaterial() {
  this->registerNeighborhood();
  this->insertQuadsInNeighborhoods(_not_ghost);
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialNonLocal<spatial_dimension>::insertQuadsInNeighborhoods(GhostType ghost_type) {

  NonLocalManager & manager = this->model->getNonLocalManager();
  UInt spatial_dimension = this->model->getSpatialDimension();
  InternalField<Real> quadrature_points_coordinates("quadrature_points_coordinates_tmp_nl", *this);
  quadrature_points_coordinates.initialize(spatial_dimension);

  /// intialize quadrature point object
  IntegrationPoint q;
  q.ghost_type = ghost_type;
  q.kind = _ek_regular;

  Mesh::type_iterator it = this->element_filter.firstType(spatial_dimension, ghost_type, _ek_regular);
  Mesh::type_iterator last_type = this->element_filter.lastType(spatial_dimension, ghost_type, _ek_regular);
  for(; it != last_type; ++it) {
    q.type = *it;
    const Array<UInt> & elem_filter = this->element_filter(*it, ghost_type);
    UInt nb_element  = elem_filter.getSize();
    if(nb_element) {
      UInt nb_quad = this->fem->getNbIntegrationPoints(*it, ghost_type);
      UInt nb_tot_quad = nb_quad * nb_element;
      
      Array<Real> & quads = quadrature_points_coordinates(*it, ghost_type);
      quads.resize(nb_tot_quad);

      this->model->getFEEngine().computeIntegrationPointsCoordinates(quads, *it, ghost_type, elem_filter);
      
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
      UInt nb_quads = this->getFEEngine().getNbIntegrationPoints(el_type, ghost_type);
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

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialNonLocal<spatial_dimension>::updateResidual(GhostType ghost_type) {
  AKANTU_EXCEPTION("this method has not been implemented");
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialNonLocal<spatial_dimension>::registerNeighborhood() {
  this->model->getNonLocalManager().registerNeighborhood(this->name, this->name);
}

/* -------------------------------------------------------------------------- */

INSTANTIATE_MATERIAL(MaterialNonLocal);

__END_AKANTU__
