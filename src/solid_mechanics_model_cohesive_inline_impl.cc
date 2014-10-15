/**
 * @file   solid_mechanics_model_cohesive_inline_impl.cc
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Fri Jan 18 09:40:01 2013
 *
 * @brief  Solid mechanics model inline implementation
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

/* -------------------------------------------------------------------------- */
inline UInt SolidMechanicsModelCohesive::getNbQuadsForFacetCheck(const Array<Element> & elements) const {
  UInt nb_quads = 0;
  UInt nb_quad_per_facet = 0;

  ElementType current_element_type = _not_defined;
  GhostType current_ghost_type = _casper;

  Array<Element>::const_iterator<Element> it  = elements.begin();
  Array<Element>::const_iterator<Element> end = elements.end();
  for (; it != end; ++it) {
    const Element & el = *it;

    if(el.type != current_element_type || el.ghost_type != current_ghost_type) {
      current_element_type = el.type;
      current_ghost_type   = el.ghost_type;

      nb_quad_per_facet = this->getFEEngine("FacetsFEEngine").getNbQuadraturePoints(el.type,
									  el.ghost_type);
    }

    nb_quads += nb_quad_per_facet;
  }

  return nb_quads;
}

/* -------------------------------------------------------------------------- */
inline UInt SolidMechanicsModelCohesive::getNbNodesPerElementList(const Array<Element> & elements) const {
  UInt nb_nodes_per_element = 0;
  UInt nb_nodes = 0;
  ElementType current_element_type = _not_defined;

  Array<Element>::const_iterator<Element> el_it = elements.begin();
  Array<Element>::const_iterator<Element> el_end = elements.end();

  for (; el_it != el_end; ++el_it) {
    const Element & el = *el_it;

    if(el.type != current_element_type) {
      current_element_type = el.type;
      nb_nodes_per_element = Mesh::getNbNodesPerElement(current_element_type);
    }

    nb_nodes += nb_nodes_per_element;
  }

  return nb_nodes;
}

/* -------------------------------------------------------------------------- */
inline UInt SolidMechanicsModelCohesive::getNbDataForElements(const Array<Element> & elements,
							      SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;

  if (elements.getSize() == 0) return size;

  /// regular element case
  if (elements(0).kind == _ek_regular) {

    switch(tag) {
    case _gst_smmc_facets: {
      size += elements.getSize() * sizeof(bool);
      break;
    }
    case _gst_smmc_facets_conn: {
      UInt nb_nodes = getNbNodesPerElementList(elements);
      size += nb_nodes * sizeof(UInt);
      break;
    }
    case _gst_smmc_facets_stress: {
      UInt nb_quads = getNbQuadsForFacetCheck(elements);
      size += nb_quads * spatial_dimension * spatial_dimension * sizeof(Real);
      break;
    }
    default: {
      size += SolidMechanicsModel::getNbDataForElements(elements, tag);
    }
    }
  }
  /// cohesive element case
  else if (elements(0).kind == _ek_cohesive) {

    switch(tag) {
    case _gst_material_id: {
      size += elements.getSize() * 2 * sizeof(UInt);
      break;
    }
    case _gst_smm_boundary: {
      UInt nb_nodes_per_element = 0;

      Array<Element>::const_iterator<Element> it  = elements.begin();
      Array<Element>::const_iterator<Element> end = elements.end();
      for (; it != end; ++it) {
	const Element & el = *it;
	nb_nodes_per_element += Mesh::getNbNodesPerElement(el.type);
      }

      // force, displacement, boundary
      size += nb_nodes_per_element * spatial_dimension * (2 * sizeof(Real) + sizeof(bool));
      break;
    }
    default: { }
    }

    if(tag != _gst_material_id && tag != _gst_smmc_facets) {
      Array<Element> * elements_per_mat = new Array<Element>[materials.size()];
      this->splitElementByMaterial(elements, elements_per_mat);

      for (UInt i = 0; i < materials.size(); ++i) {
	size += materials[i]->getNbDataForElements(elements_per_mat[i], tag);
      }
      delete [] elements_per_mat;
    }
  }

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
inline void SolidMechanicsModelCohesive::packElementData(CommunicationBuffer & buffer,
							 const Array<Element> & elements,
							 SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  if (elements.getSize() == 0) return;

  if (elements(0).kind == _ek_regular) {

    switch(tag) {

    case _gst_smmc_facets: {
      packElementalDataHelper(inserter->getInsertionFacetsByElement(),
			      buffer, elements, false, getFEEngine());
      break;
    }
    case _gst_smmc_facets_conn: {
      packElementalDataHelper(*global_connectivity, buffer, elements,
			      false, getFEEngine());
      break;
    }
    case _gst_smmc_facets_stress: {
      packFacetStressDataHelper(facet_stress, buffer, elements);
      break;
    }
    default: {
      SolidMechanicsModel::packElementData(buffer, elements, tag);
    }
    }
  }
  else if (elements(0).kind == _ek_cohesive) {

    switch(tag) {

    case _gst_material_id: {
      packElementalDataHelper(element_index_by_material, buffer,
			      elements, false, getFEEngine("CohesiveFEEngine"));
      break;
    }
    case _gst_smm_boundary: {
      packNodalDataHelper(*force, buffer, elements, mesh);
      packNodalDataHelper(*velocity, buffer, elements, mesh);
      packNodalDataHelper(*blocked_dofs, buffer, elements, mesh);
      break;
    }
    default: { }
    }

    if(tag != _gst_material_id && tag != _gst_smmc_facets) {
      Array<Element> * elements_per_mat = new Array<Element>[materials.size()];
      splitElementByMaterial(elements, elements_per_mat);

      for (UInt i = 0; i < materials.size(); ++i) {
	materials[i]->packElementData(buffer, elements_per_mat[i], tag);
      }

      delete [] elements_per_mat;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void SolidMechanicsModelCohesive::unpackElementData(CommunicationBuffer & buffer,
							   const Array<Element> & elements,
							   SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  if (elements.getSize() == 0) return;

  if (elements(0).kind == _ek_regular) {

    switch(tag) {
    case _gst_smmc_facets: {
      unpackElementalDataHelper(inserter->getInsertionFacetsByElement(),
				buffer, elements, false, getFEEngine());
      break;
    }
    case _gst_smmc_facets_conn: {
      unpackElementalDataHelper(*global_connectivity, buffer, elements,
				false, getFEEngine());
      break;
    }
    case _gst_smmc_facets_stress: {
      unpackFacetStressDataHelper(facet_stress, buffer, elements);
      break;
    }
    default: {
      SolidMechanicsModel::unpackElementData(buffer, elements, tag);
    }
    }
  }
  else if (elements(0).kind == _ek_cohesive) {

    switch(tag) {
    case _gst_material_id: {
      unpackElementalDataHelper(element_index_by_material, buffer,
				elements, false, getFEEngine("CohesiveFEEngine"));
      break;
    }
    case _gst_smm_boundary: {
      unpackNodalDataHelper(*force, buffer, elements, mesh);
      unpackNodalDataHelper(*velocity, buffer, elements, mesh);
      unpackNodalDataHelper(*blocked_dofs, buffer, elements, mesh);
      break;
    }
    default: { }
    }

    if(tag != _gst_material_id && tag != _gst_smmc_facets) {
      Array<Element> * elements_per_mat = new Array<Element>[materials.size()];
      splitElementByMaterial(elements, elements_per_mat);

      for (UInt i = 0; i < materials.size(); ++i) {
	materials[i]->unpackElementData(buffer, elements_per_mat[i], tag);
      }

      delete [] elements_per_mat;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void SolidMechanicsModelCohesive::packFacetStressDataHelper(const ElementTypeMapArray<T> & data_to_pack,
								   CommunicationBuffer & buffer,
								   const Array<Element> & elements) const {
  packUnpackFacetStressDataHelper<T, true>(const_cast<ElementTypeMapArray<T> &>(data_to_pack),
					   buffer, elements);
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void SolidMechanicsModelCohesive::unpackFacetStressDataHelper(ElementTypeMapArray<T> & data_to_unpack,
								     CommunicationBuffer & buffer,
								     const Array<Element> & elements) const {
  packUnpackFacetStressDataHelper<T, false>(data_to_unpack, buffer, elements);
}

/* -------------------------------------------------------------------------- */
template<typename T, bool pack_helper>
inline void SolidMechanicsModelCohesive::packUnpackFacetStressDataHelper(ElementTypeMapArray<T> & data_to_pack,
									 CommunicationBuffer & buffer,
									 const Array<Element> & element) const {
  ElementType current_element_type = _not_defined;
  GhostType current_ghost_type = _casper;
  UInt nb_quad_per_elem = 0;
  UInt sp2 = spatial_dimension * spatial_dimension;
  UInt nb_component = sp2 * 2;
  bool element_rank = 0;
  Mesh & mesh_facets = inserter->getMeshFacets();

  Array<T> * vect = NULL;
  Array<std::vector<Element> > * element_to_facet = NULL;

  Array<Element>::const_iterator<Element> it  = element.begin();
  Array<Element>::const_iterator<Element> end = element.end();
  for (; it != end; ++it) {
    const Element & el = *it;
    if(el.type != current_element_type || el.ghost_type != current_ghost_type) {
      current_element_type = el.type;
      current_ghost_type   = el.ghost_type;
      vect = &data_to_pack(el.type, el.ghost_type);

      element_to_facet =
	&( mesh_facets.getElementToSubelement(el.type, el.ghost_type) );

      nb_quad_per_elem = this->getFEEngine("FacetsFEEngine").getNbQuadraturePoints(el.type,
									 el.ghost_type);
    }

    if (pack_helper)
      element_rank = (*element_to_facet)(el.element)[0].ghost_type != _not_ghost;
    else
      element_rank = (*element_to_facet)(el.element)[0].ghost_type == _not_ghost;

    for (UInt q = 0; q < nb_quad_per_elem; ++q) {
      Vector<T> data(vect->storage()
		     + (el.element * nb_quad_per_elem + q) * nb_component
		     + element_rank * sp2,
		     sp2);

      if(pack_helper)
	buffer << data;
      else
	buffer >> data;
    }
  }
}
