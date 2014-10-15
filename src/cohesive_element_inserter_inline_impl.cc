/**
 * @file   cohesive_element_inserter_inline_impl.cc
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Tue Dec  3 15:20:53 2013
 *
 * @brief  Cohesive element inserter inline functions
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
#include "cohesive_element_inserter.hh"

/* -------------------------------------------------------------------------- */
inline UInt CohesiveElementInserter::getNbDataForElements(const Array<Element> & elements,
							  SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;

  if (tag == _gst_ce_inserter) {
    UInt nb_nodes = 0;

    Array<Element>::const_iterator<Element> it  = elements.begin();
    Array<Element>::const_iterator<Element> end = elements.end();
    for (; it != end; ++it) {
      const Element & el = *it;
      nb_nodes += Mesh::getNbNodesPerElement(el.type);
    }

    size += nb_nodes * sizeof(UInt);
  }

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
inline void CohesiveElementInserter::packElementData(CommunicationBuffer & buffer,
						     const Array<Element> & elements,
						     SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  if (tag == _gst_ce_inserter)
    packUnpackGlobalConnectivity<true>(buffer, elements);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void CohesiveElementInserter::unpackElementData(CommunicationBuffer & buffer,
						       const Array<Element> & elements,
						       SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  if (tag == _gst_ce_inserter)
    packUnpackGlobalConnectivity<false>(buffer, elements);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<bool pack_mode>
inline void CohesiveElementInserter::packUnpackGlobalConnectivity(CommunicationBuffer & buffer,
								  const Array<Element> & elements) const {
  AKANTU_DEBUG_IN();

  ElementType current_element_type = _not_defined;
  GhostType current_ghost_type = _casper;

  Array<UInt>::iterator<Vector<UInt> > conn_begin;
  UInt nb_nodes_per_elem = 0;
  UInt index;

  Array<Element>::const_iterator<Element> it  = elements.begin();
  Array<Element>::const_iterator<Element> end = elements.end();
  for (; it != end; ++it) {
    const Element & el = *it;

    if (el.type != current_element_type || el.ghost_type != current_ghost_type) {
      current_element_type = el.type;
      current_ghost_type   = el.ghost_type;

      nb_nodes_per_elem = Mesh::getNbNodesPerElement(current_element_type);

      conn_begin = mesh.connectivities(current_element_type,
				       current_ghost_type).begin(nb_nodes_per_elem);
    }

    /// get element connectivity
    Vector<UInt> & current_conn = conn_begin[el.element];

    /// loop on all connectivity nodes
    for (UInt n = 0; n < nb_nodes_per_elem; ++n) {
      UInt node = current_conn(n);

      if (pack_mode) {
	/// if node is local or master pack its global id, otherwise
	/// dummy data
	if (mesh.isLocalOrMasterNode(node))
	  index = mesh.getNodeGlobalId(node);
	else
	  index = UInt(-1);

	buffer << index;
      }
      else {
	buffer >> index;

	/// update slave nodes' index
	if (index != UInt(-1) && mesh.isSlaveNode(node))
	  (*mesh.nodes_global_ids)(node) = index;
      }
    }
  }

  AKANTU_DEBUG_OUT();
}
