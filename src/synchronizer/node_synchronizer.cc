/**
 * @file   node_synchronizer.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Nov 15 2017
 *
 * @brief  Implementation of the node synchronizer
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
#include "node_synchronizer.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
NodeSynchronizer::NodeSynchronizer(Mesh & mesh, const ID & id,
                                   const bool register_to_event_manager,
                                   EventHandlerPriority event_priority)
    : SynchronizerImpl<UInt>(mesh.getCommunicator(), id),
      mesh(mesh) {
  AKANTU_DEBUG_IN();

  if (register_to_event_manager) {
    this->mesh.registerEventHandler(*this, event_priority);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
NodeSynchronizer::~NodeSynchronizer() = default;

/* -------------------------------------------------------------------------- */
Int NodeSynchronizer::getRank(const UInt & node) const {
  return this->mesh.getNodePrank(node);
}

/* -------------------------------------------------------------------------- */
void NodeSynchronizer::onNodesAdded(const Array<UInt> & /*nodes_list*/,
                                    const NewNodesEvent & /*unused*/) {
  std::map<UInt, std::vector<UInt>> nodes_per_proc;

  // recreates fully the schemes due to changes of global ids
  // \TODO add an event to handle global id changes
  for (auto && data : communications.iterateSchemes(_recv)) {
    auto & scheme = data.second;
    scheme.resize(0);
  }

  for (auto && local_id : arange(mesh.getNbNodes())) {
    if (not mesh.isSlaveNode(local_id)) {
      continue; // local, master or pure ghost
    }

    auto global_id = mesh.getNodeGlobalId(local_id);
    auto proc = mesh.getNodePrank(local_id);
    AKANTU_DEBUG_ASSERT(
        proc != -1,
        "The node " << local_id << " does not have a valid associated prank");
    nodes_per_proc[proc].push_back(global_id);

    auto & scheme = communications.createScheme(proc, _recv);
    scheme.push_back(local_id);
  }

  std::vector<CommunicationRequest> send_requests;
  for (auto && pair : communications.iterateSchemes(_recv)) {
    auto proc = pair.first;
    AKANTU_DEBUG_ASSERT(proc != UInt(-1),
                        "For real I should send something to proc -1");

    // if proc not in nodes_per_proc this should insert an empty array to send
    send_requests.push_back(communicator.asyncSend(
        nodes_per_proc[proc], proc, Tag::genTag(rank, proc, 0xcafe)));
  }

  for (auto && data : communications.iterateSchemes(_send)) {
    auto proc = data.first;
    auto & scheme = data.second;
    CommunicationStatus status;

    auto tag = Tag::genTag(proc, rank, 0xcafe);
    communicator.probe<UInt>(proc, tag, status);

    scheme.resize(status.size());
    communicator.receive(scheme, proc, tag);
    std::transform(scheme.begin(), scheme.end(), scheme.begin(),
                   [&](auto & gnode) { return mesh.getNodeLocalId(gnode); });
  }

  // communicator.receiveAnyNumber<UInt>(
  //     send_requests,
  //     [&](auto && proc, auto && nodes) {
  //       auto & scheme = communications.createScheme(proc, _send);
  //       scheme.resize(nodes.size());
  //       for (auto && data : enumerate(nodes)) {
  //         auto global_id = std::get<1>(data);
  //         auto local_id = mesh.getNodeLocalId(global_id);
  //         AKANTU_DEBUG_ASSERT(local_id != UInt(-1),
  //                             "The global node " << global_id
  //                                                << "is not known on rank "
  //                                                << rank);
  //         scheme[std::get<0>(data)] = local_id;
  //       }
  //     },
  //     Tag::genTag(rank, count, 0xcafe));
  // ++count;

  Communicator::waitAll(send_requests);
  Communicator::freeCommunicationRequest(send_requests);

  this->entities_changed = true;
}

/* -------------------------------------------------------------------------- */
UInt NodeSynchronizer::sanityCheckDataSize(const Array<UInt> & nodes,
                                           const SynchronizationTag & tag,
                                           bool from_comm_desc) const {
  UInt size =
      SynchronizerImpl<UInt>::sanityCheckDataSize(nodes, tag, from_comm_desc);

  // global id
  if (tag != SynchronizationTag::_giu_global_conn) {
    size += sizeof(UInt) * nodes.size();
  }

  // flag
  size += sizeof(NodeFlag) * nodes.size();

  // positions
  size += mesh.getSpatialDimension() * sizeof(Real) * nodes.size();

  return size;
}

/* -------------------------------------------------------------------------- */
void NodeSynchronizer::packSanityCheckData(
    CommunicationBuffer & buffer, const Array<UInt> & nodes,
    const SynchronizationTag & tag) const {
  auto dim = mesh.getSpatialDimension();
  for (auto && node : nodes) {
    if (tag != SynchronizationTag::_giu_global_conn) {
      buffer << mesh.getNodeGlobalId(node);
    }
    buffer << mesh.getNodeFlag(node);
    buffer << Vector<Real>(mesh.getNodes().begin(dim)[node]);
  }
}

/* -------------------------------------------------------------------------- */
void NodeSynchronizer::unpackSanityCheckData(CommunicationBuffer & buffer,
                                             const Array<UInt> & nodes,
                                             const SynchronizationTag & tag,
                                             UInt proc, UInt rank) const {
  auto dim = mesh.getSpatialDimension();

#ifndef AKANTU_NDEBUG
  auto periodic = [&](auto && flag) { return flag & NodeFlag::_periodic_mask; };
  auto distrib = [&](auto && flag) { return flag & NodeFlag::_shared_mask; };
#endif

  for (auto && node : nodes) {
    if (tag != SynchronizationTag::_giu_global_conn) {
      UInt global_id;
      buffer >> global_id;
      AKANTU_DEBUG_ASSERT(global_id == mesh.getNodeGlobalId(node),
                          "The nodes global ids do not match: "
                              << global_id
                              << " != " << mesh.getNodeGlobalId(node));
    }

    NodeFlag flag;
    buffer >> flag;
    AKANTU_DEBUG_ASSERT(
        (periodic(flag) == periodic(mesh.getNodeFlag(node))) and
            (((distrib(flag) == NodeFlag::_master) and
              (distrib(mesh.getNodeFlag(node)) ==
               NodeFlag::_slave)) or // master to slave
             ((distrib(flag) == NodeFlag::_slave) and
              (distrib(mesh.getNodeFlag(node)) ==
               NodeFlag::_master)) or // reverse comm slave to master
             (distrib(mesh.getNodeFlag(node)) ==
                  NodeFlag::_pure_ghost or // pure ghost nodes
              distrib(flag) == NodeFlag::_pure_ghost)),
        "The node flags: " << flag << " and " << mesh.getNodeFlag(node));

    Vector<Real> pos_remote(dim);
    buffer >> pos_remote;
    Vector<Real> pos(mesh.getNodes().begin(dim)[node]);

    auto dist = pos_remote.distance(pos);
    if (not Math::are_float_equal(dist, 0.)) {
      AKANTU_EXCEPTION("Unpacking an unknown value for the node "
                       << node << "(position " << pos << " != buffer "
                       << pos_remote << ") [" << dist << "] - tag: " << tag
                       << " comm from " << proc << " to " << rank);
    }
  }
}

/* -------------------------------------------------------------------------- */
void NodeSynchronizer::fillEntityToSend(Array<UInt> & nodes_to_send) {
  UInt nb_nodes = mesh.getNbNodes();

  this->entities_from_root.clear();
  nodes_to_send.resize(0);

  for (UInt n : arange(nb_nodes)) {
    if (not mesh.isLocalOrMasterNode(n)) {
      continue;
    }

    entities_from_root.push_back(n);
  }

  for (auto n : entities_from_root) {
    UInt global_node = mesh.getNodeGlobalId(n);
    nodes_to_send.push_back(global_node);
  }
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
