/**
 * @file   node_group.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Thu Feb 01 2018
 *
 * @brief  Implementation of the node group
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
#include "node_group.hh"
#include "dumpable.hh"
#include "dumpable_inline_impl.hh"
#include "mesh.hh"
#if defined(AKANTU_USE_IOHELPER)
#include "dumper_iohelper_paraview.hh"
#endif
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
NodeGroup::NodeGroup(const std::string & name, const Mesh & mesh,
                     const std::string & id, const MemoryID & memory_id)
    : Memory(id, memory_id), name(name),
      node_group(alloc<UInt>(std::string(this->id + ":nodes"), 0, 1)) {
#if defined(AKANTU_USE_IOHELPER)
  this->registerDumper<DumperParaview>("paraview_" + name, name, true);
  auto field = std::make_shared<dumpers::NodalField<Real, true>>(
      mesh.getNodes(), 0, 0, &this->getNodes());
  this->getDumper().registerField("positions", field);
#endif
}

/* -------------------------------------------------------------------------- */
NodeGroup::~NodeGroup() = default;

/* -------------------------------------------------------------------------- */
void NodeGroup::clear() { node_group.resize(0); }
/* -------------------------------------------------------------------------- */
//bool NodeGroup::empty() { return node_group.empty(); }
/* -------------------------------------------------------------------------- */
void NodeGroup::optimize() {
  std::sort(node_group.begin(), node_group.end());
  Array<UInt>::iterator<> end =
      std::unique(node_group.begin(), node_group.end());
  node_group.resize(end - node_group.begin());
}

/* -------------------------------------------------------------------------- */
void NodeGroup::append(const NodeGroup & other_group) {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = node_group.size();

  /// append new nodes to current list
  node_group.resize(nb_nodes + other_group.node_group.size());
  std::copy(other_group.node_group.begin(), other_group.node_group.end(),
            node_group.begin() + nb_nodes);

  optimize();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NodeGroup::printself(std::ostream & stream, int indent) const {
  std::string space(indent, AKANTU_INDENT);

  stream << space << "NodeGroup [" << std::endl;
  stream << space << " + name: " << name << std::endl;
  node_group.printself(stream, indent + 1);
  stream << space << "]" << std::endl;
}

} // namespace akantu
