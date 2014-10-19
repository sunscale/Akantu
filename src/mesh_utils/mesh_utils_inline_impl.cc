/**
 * @file   mesh_utils_inline_impl.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Aug 20 2010
 * @date last modification: Mon Jun 09 2014
 *
 * @brief  Mesh utils inline functions
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
inline bool MeshUtils::hasElement(const Array<UInt> & connectivity,
				  const Element & el,
				  const Vector<UInt> & nodes) {

  UInt nb_nodes_per_element = connectivity.getNbComponent();

  const Vector<UInt> el_nodes(connectivity.storage()
			      + el.element * nb_nodes_per_element,
			      nb_nodes_per_element);

  UInt matched_nodes = 0;

  for (UInt n = 0; n < nodes.size(); ++n) {
    if (std::find(el_nodes.storage(),
		  el_nodes.storage() + nb_nodes_per_element,
		  nodes[n])
	!= (el_nodes.storage() + nb_nodes_per_element))
      ++matched_nodes;
  }

  return (matched_nodes == nodes.size());
}

/* -------------------------------------------------------------------------- */
inline void MeshUtils::updateElementalConnectivity(Mesh & mesh,
						   UInt old_node,
						   UInt new_node,
						   const std::vector<Element> & element_list,
						   const std::vector<Element> * facet_list) {
  AKANTU_DEBUG_IN();

  ElementType el_type = _not_defined;
  GhostType gt_type = _casper;
  Array<UInt> * conn_elem = NULL;
#if defined(AKANTU_COHESIVE_ELEMENT)
  const Array<Element> * cohesive_facets = NULL;
#endif
  UInt nb_nodes_per_element = 0;
  UInt * n_update = NULL;

  for (UInt el = 0; el < element_list.size(); ++el) {
    const Element & elem = element_list[el];
    if (elem.type == _not_defined) continue;

    if (elem.type != el_type || elem.ghost_type != gt_type) {
      el_type = elem.type;
      gt_type = elem.ghost_type;
      conn_elem = & mesh.getConnectivity(el_type, gt_type);
      nb_nodes_per_element = conn_elem->getNbComponent();
#if defined(AKANTU_COHESIVE_ELEMENT)
      if (elem.kind == _ek_cohesive)
	cohesive_facets = & mesh.getMeshFacets().getSubelementToElement(el_type, gt_type);
#endif
    }

#if defined(AKANTU_COHESIVE_ELEMENT)
    if (elem.kind == _ek_cohesive) {

      AKANTU_DEBUG_ASSERT(facet_list != NULL,
			  "Provide a facet list in order to update cohesive elements");

      /// loop over cohesive element's facets
      for (UInt f = 0, n = 0; f < 2; ++f, n += nb_nodes_per_element / 2) {
	const Element & facet = (*cohesive_facets)(elem.element, f);

	/// skip facets if not present in the list
	if (std::find(facet_list->begin(), facet_list->end(), facet)
	    == facet_list->end()) continue;

	n_update
	  = std::find(conn_elem->storage() + elem.element * nb_nodes_per_element + n,
		      conn_elem->storage() + elem.element * nb_nodes_per_element + n
		      + nb_nodes_per_element / 2,
		      old_node);

	AKANTU_DEBUG_ASSERT(n_update != conn_elem->storage()
			    + elem.element * nb_nodes_per_element + n
			    + nb_nodes_per_element / 2,
			    "Node not found in current element");

	/// update connectivity
	*n_update = new_node;
      }
    }
    else {
#endif
      n_update
	= std::find(conn_elem->storage() + elem.element * nb_nodes_per_element,
		    conn_elem->storage() + elem.element * nb_nodes_per_element
		    + nb_nodes_per_element,
		    old_node);

      AKANTU_DEBUG_ASSERT(n_update != conn_elem->storage()
			  + elem.element * nb_nodes_per_element
			  + nb_nodes_per_element,
			  "Node not found in current element");

      /// update connectivity
      *n_update = new_node;
#if defined(AKANTU_COHESIVE_ELEMENT)
    }
#endif
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline bool MeshUtils::removeElementsInVector(const std::vector<Element> & elem_to_remove,
					      std::vector<Element> & elem_list) {
  if (elem_list.size() <= elem_to_remove.size())
    return true;

  std::vector<Element>::const_iterator el_it = elem_to_remove.begin();
  std::vector<Element>::const_iterator el_last = elem_to_remove.end();
  std::vector<Element>::iterator el_del;

  UInt deletions = 0;

  for (; el_it != el_last; ++el_it) {
    el_del = std::find(elem_list.begin(), elem_list.end(), *el_it);

    if (el_del != elem_list.end()) {
      elem_list.erase(el_del);
      ++deletions;
    }
  }

  AKANTU_DEBUG_ASSERT(deletions == 0 || deletions == elem_to_remove.size(),
		      "Not all elements have been erased");

  return deletions == 0;
}
