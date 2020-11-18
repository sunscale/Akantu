/**
 * @file   node_group_inline_impl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Sun Aug 13 2017
 *
 * @brief  Node group inline function definitions
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

namespace akantu {

/* -------------------------------------------------------------------------- */
inline NodeGroup::const_node_iterator NodeGroup::begin() const {
  return node_group.begin();
}

/* -------------------------------------------------------------------------- */
inline NodeGroup::const_node_iterator NodeGroup::end() const {
  return node_group.end();
}

/* -------------------------------------------------------------------------- */
inline NodeGroup::const_node_iterator NodeGroup::add(UInt node,
                                                     bool check_for_duplicate) {
  if (check_for_duplicate) {
    const_node_iterator it = std::find(begin(), end(), node);
    if (it == node_group.end()) {
      node_group.push_back(node);
      return (node_group.end() - 1);
    }
    return it;
  }

  node_group.push_back(node);
  return (node_group.end() - 1);
}

/* -------------------------------------------------------------------------- */
inline void NodeGroup::remove(UInt node) {
  Array<UInt>::iterator<> it = this->node_group.begin();
  Array<UInt>::iterator<> end = this->node_group.end();
  AKANTU_DEBUG_ASSERT(it != end, "The node group is empty!!");
  for (; it != node_group.end(); ++it) {
    if (*it == node) {
      it = node_group.erase(it);
    }
  }
  AKANTU_DEBUG_ASSERT(it != end, "The node was not found!");
}

/* -------------------------------------------------------------------------- */
inline UInt NodeGroup::size() const { return node_group.size(); }

/* -------------------------------------------------------------------------- */
struct FilterFunctor;

template <typename T> void NodeGroup::applyNodeFilter(T & filter) {
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

} // namespace akantu
