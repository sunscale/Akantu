/**
 * @file   element_group_inline_impl.hh
 *
 * @author Dana Christen <dana.christen@gmail.com>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Nov 13 2013
 * @date last modification: Sun Aug 13 2017
 *
 * @brief  Stores information relevent to the notion of domain boundary and
 * surfaces.
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "element_group.hh"
#include "mesh.hh"

/* -------------------------------------------------------------------------- */

#ifndef AKANTU_ELEMENT_GROUP_INLINE_IMPL_HH_
#define AKANTU_ELEMENT_GROUP_INLINE_IMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
inline void ElementGroup::add(const Element & el, bool add_nodes,
                              bool check_for_duplicate) {
  this->add(el.type, el.element, el.ghost_type, add_nodes, check_for_duplicate);
}

/* -------------------------------------------------------------------------- */

inline void ElementGroup::add(ElementType type, UInt element,
                              GhostType ghost_type, bool add_nodes,
                              bool check_for_duplicate) {

  addElement(type, element, ghost_type);

  if (add_nodes) {
    Array<UInt>::const_vector_iterator it =
        mesh.getConnectivity(type, ghost_type)
            .begin(mesh.getNbNodesPerElement(type)) +
        element;
    const Vector<UInt> & conn = *it;
    for (UInt i = 0; i < conn.size(); ++i) {
      addNode(conn[i], check_for_duplicate);
    }
  }
}

/* -------------------------------------------------------------------------- */
inline void ElementGroup::addNode(UInt node_id, bool check_for_duplicate) {
  node_group.add(node_id, check_for_duplicate);
}

/* -------------------------------------------------------------------------- */
inline void ElementGroup::removeNode(UInt node_id) {
  node_group.remove(node_id);
}

/* -------------------------------------------------------------------------- */
inline void ElementGroup::addElement(ElementType elem_type,
                                     UInt elem_id,
                                     GhostType ghost_type) {
  if (!(elements.exists(elem_type, ghost_type))) {
    elements.alloc(0, 1, elem_type, ghost_type);
  }

  elements(elem_type, ghost_type).push_back(elem_id);
  this->dimension = UInt(
      std::max(Int(this->dimension), Int(mesh.getSpatialDimension(elem_type))));
}

/* -------------------------------------------------------------------------- */
inline UInt ElementGroup::getNbNodes() const { return node_group.size(); }

/* -------------------------------------------------------------------------- */
inline ElementGroup::type_iterator
ElementGroup::firstType(UInt dim, GhostType ghost_type,
                        ElementKind kind) const {
  return elements.elementTypes(dim, ghost_type, kind).begin();
}

/* -------------------------------------------------------------------------- */
inline ElementGroup::type_iterator
ElementGroup::lastType(UInt dim, GhostType ghost_type,
                       ElementKind kind) const {
  return elements.elementTypes(dim, ghost_type, kind).end();
}

/* -------------------------------------------------------------------------- */
inline ElementGroup::const_element_iterator
ElementGroup::begin(ElementType type,
                    GhostType ghost_type) const {
  if (elements.exists(type, ghost_type)) {
    return elements(type, ghost_type).begin();
  }
  return empty_elements.begin();
}

/* -------------------------------------------------------------------------- */
inline ElementGroup::const_element_iterator
ElementGroup::end(ElementType type,
                  GhostType ghost_type) const {
  if (elements.exists(type, ghost_type)) {
    return elements(type, ghost_type).end();
  }
  return empty_elements.end();
}

/* -------------------------------------------------------------------------- */
inline const Array<UInt> &
ElementGroup::getElements(ElementType type,
                          GhostType ghost_type) const {
  if (elements.exists(type, ghost_type)) {
    return elements(type, ghost_type);
  }
  return empty_elements;
}

/* -------------------------------------------------------------------------- */

} // namespace akantu

#endif /* AKANTU_ELEMENT_GROUP_INLINE_IMPL_HH_ */
