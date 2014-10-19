/**
 * @file   node_group_inline_impl.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Nov 13 2013
 * @date last modification: Tue Sep 02 2014
 *
 * @brief  Node group inline function definitions
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

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
inline NodeGroup::const_node_iterator NodeGroup::begin() const {
  return node_group.begin();
}

/* -------------------------------------------------------------------------- */
inline NodeGroup::const_node_iterator NodeGroup::end() const {
  return node_group.end();
}

/* -------------------------------------------------------------------------- */
inline NodeGroup::const_node_iterator NodeGroup::add(UInt node, bool check_for_duplicate) {
  if(check_for_duplicate) {
    const_node_iterator it = std::find(begin(), end(), node);
    if(it == node_group.end()) {
      node_group.push_back(node);
      return (node_group.end() - 1);
    }
    return it;
  } else {
    node_group.push_back(node);
    return (node_group.end() - 1);
  }
}

/* -------------------------------------------------------------------------- */
inline UInt NodeGroup::getSize() const {
  return node_group.getSize();
}

/* -------------------------------------------------------------------------- */
class FilterFunctor;

template <typename T>
void NodeGroup::applyNodeFilter(T & filter) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(T::type == FilterFunctor::_node_filter_functor,
		      "NodeFilter can only apply node filter functor");

  Array<UInt>::iterator<> it = this->node_group.begin();

  for (; it != node_group.end(); ++it) {
    /// filter == true -> keep node
    if (!filter(*it)) { 
      it = node_group.erase(it);
    }
  }

  AKANTU_DEBUG_OUT();
}


__END_AKANTU__
