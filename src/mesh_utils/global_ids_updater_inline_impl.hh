/**
 * @file   global_ids_updater_inline_impl.hh
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Oct 02 2015
 * @date last modification: Sun Aug 13 2017
 *
 * @brief  Implementation of the inline functions of GlobalIdsUpdater
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "communicator.hh"
#include "global_ids_updater.hh"
#include "mesh.hh"
#include "mesh_accessor.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_GLOBAL_IDS_UPDATER_INLINE_IMPL_HH_
#define AKANTU_GLOBAL_IDS_UPDATER_INLINE_IMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
inline UInt GlobalIdsUpdater::getNbData(const Array<Element> & elements,
                                        const SynchronizationTag & tag) const {
  UInt size = 0;
  if (tag == SynchronizationTag::_giu_global_conn) {
    size +=
        Mesh::getNbNodesPerElementList(elements) * sizeof(UInt) + sizeof(int);
#ifndef AKANTU_NDEBUG
    size += sizeof(NodeFlag) * Mesh::getNbNodesPerElementList(elements);
#endif
  }
  return size;
}

/* -------------------------------------------------------------------------- */
inline void GlobalIdsUpdater::packData(CommunicationBuffer & buffer,
                                       const Array<Element> & elements,
                                       const SynchronizationTag & tag) const {
  if (tag != SynchronizationTag::_giu_global_conn) {
    return;
  }

  int prank = mesh.getCommunicator().whoAmI();

  auto & global_nodes_ids = mesh.getGlobalNodesIds();
  buffer << prank;

  for (const auto & element : elements) {
    /// get element connectivity
    const Vector<UInt> current_conn =
        const_cast<const Mesh &>(mesh).getConnectivity(element);

    /// loop on all connectivity nodes
    for (auto node : current_conn) {
      UInt index = -1;
      if ((this->reduce and mesh.isLocalOrMasterNode(node)) or
          (not this->reduce and not mesh.isPureGhostNode(node))) {
        index = global_nodes_ids(node);
      }
      buffer << index;
#ifndef AKANTU_NDEBUG
      buffer << mesh.getNodeFlag(node);
#endif
    }
  }
}

/* -------------------------------------------------------------------------- */
inline void GlobalIdsUpdater::unpackData(CommunicationBuffer & buffer,
                                         const Array<Element> & elements,
                                         const SynchronizationTag & tag) {
  if (tag != SynchronizationTag::_giu_global_conn) {
    return;
  }

  MeshAccessor mesh_accessor(mesh);
  auto & global_nodes_ids = mesh_accessor.getNodesGlobalIds();

  int proc;
  buffer >> proc;

  for (const auto & element : elements) {
    /// get element connectivity
    Vector<UInt> current_conn =
        const_cast<const Mesh &>(mesh).getConnectivity(element);

    /// loop on all connectivity nodes
    for (auto node : current_conn) {
      UInt index;
      buffer >> index;
#ifndef AKANTU_NDEBUG
      NodeFlag node_flag;
      buffer >> node_flag;
      if (reduce) {
        nodes_flags[node].push_back(std::make_pair(proc, node_flag));
      }
#endif

      if (index == UInt(-1)) {
        continue;
      }

      if (mesh.isSlaveNode(node)) {
        auto & gid = global_nodes_ids(node);
        AKANTU_DEBUG_ASSERT(gid == UInt(-1) or gid == index,
                            "The node already has a global id, from proc "
                                << proc << ", different from the one received "
                                << gid << " " << index);
        gid = index;
        mesh_accessor.setNodePrank(node, proc);
      }
    }
  }
}

} // namespace akantu

#endif /* AKANTU_GLOBAL_IDS_UPDATER_INLINE_IMPL_HH_ */
