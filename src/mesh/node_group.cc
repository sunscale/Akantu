/**
 * @file   node_group.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Nov 13 2013
 * @date last modification: Mon Jun 09 2014
 *
 * @brief  Implementation of the node group
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

#include "node_group.hh"


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
NodeGroup::NodeGroup(const std::string & name,
                     const std::string & id,
                     const MemoryID & memory_id) :
  Memory(id, memory_id),
  name(name),
  node_group(alloc<UInt>(id + ":nodes", 0, 1)) {
}

/* -------------------------------------------------------------------------- */
NodeGroup::~NodeGroup() {}

/* -------------------------------------------------------------------------- */
void NodeGroup::empty() {
  node_group.resize(0);
}

/* -------------------------------------------------------------------------- */
void NodeGroup::optimize() {
  std::sort(node_group.begin(), node_group.end());
  Array<UInt>::iterator<> end = std::unique(node_group.begin(), node_group.end());
  node_group.resize(end - node_group.begin());
}

/* -------------------------------------------------------------------------- */
void NodeGroup::append(const NodeGroup & other_group) {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = node_group.getSize();

  /// append new nodes to current list
  node_group.resize(nb_nodes + other_group.node_group.getSize());
  std::copy(other_group.node_group.begin(),
	    other_group.node_group.end(),
	    node_group.begin() + nb_nodes);

  optimize();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NodeGroup::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "NodeGroup [" << std::endl;
  stream << space << " + name: " << name << std::endl;
  node_group.printself(stream, indent + 1);
  stream << space << "]" << std::endl;
}


__END_AKANTU__
