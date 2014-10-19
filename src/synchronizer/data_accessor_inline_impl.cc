/**
 * @file   data_accessor_inline_impl.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Sep 18 2013
 * @date last modification: Thu Jun 05 2014
 *
 * @brief  Implementation of the inline functions of the DataAccessor class
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
template<typename T, bool pack_helper>
inline void DataAccessor::packUnpackNodalDataHelper(Array<T> & data,
						    CommunicationBuffer & buffer,
						    const Array<Element> & elements,
						    const Mesh & mesh) {
  UInt nb_component = data.getNbComponent();
  UInt nb_nodes_per_element = 0;

  ElementType current_element_type = _not_defined;
  GhostType current_ghost_type = _casper;
  UInt * conn = NULL;

  Array<Element>::const_iterator<Element> it  = elements.begin();
  Array<Element>::const_iterator<Element> end = elements.end();
  for (; it != end; ++it) {
    const Element & el = *it;
    if(el.type != current_element_type || el.ghost_type != current_ghost_type) {
      current_element_type = el.type;
      current_ghost_type   = el.ghost_type;
      conn = mesh.getConnectivity(el.type, el.ghost_type).storage();
      nb_nodes_per_element = Mesh::getNbNodesPerElement(el.type);
    }

    UInt el_offset  = el.element * nb_nodes_per_element;
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      UInt offset_conn = conn[el_offset + n];
      Vector<T> data_vect(data.storage() + offset_conn * nb_component,
			  nb_component);

      if(pack_helper)
	buffer << data_vect;
      else
	buffer >> data_vect;
    }
  }
}


/* -------------------------------------------------------------------------- */
template<typename T, bool pack_helper>
inline void DataAccessor::packUnpackElementalDataHelper(ElementTypeMapArray<T> & data_to_pack,
							CommunicationBuffer & buffer,
							const Array<Element> & element,
							bool per_quadrature_point_data,
							const FEEngine & fem) {
  ElementType current_element_type = _not_defined;
  GhostType current_ghost_type = _casper;
  UInt nb_quad_per_elem = 0;
  UInt nb_component = 0;

  Array<T> * vect = NULL;

  Array<Element>::const_iterator<Element> it  = element.begin();
  Array<Element>::const_iterator<Element> end = element.end();
  for (; it != end; ++it) {
    const Element & el = *it;
    if(el.type != current_element_type || el.ghost_type != current_ghost_type) {
      current_element_type = el.type;
      current_ghost_type   = el.ghost_type;
      vect = &data_to_pack(el.type, el.ghost_type);
      if(per_quadrature_point_data)
        nb_quad_per_elem = fem.getNbQuadraturePoints(el.type,
						     el.ghost_type);
      else nb_quad_per_elem = 1;
      nb_component = vect->getNbComponent();
    }

    Vector<T> data(vect->storage() + el.element * nb_component * nb_quad_per_elem,
		   nb_component * nb_quad_per_elem);
    if(pack_helper)
      buffer << data;
    else
      buffer >> data;
  }
}
