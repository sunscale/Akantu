/**
 * @file   node_info_per_processor.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Mar 16 2016
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  Please type the brief for file: Helper classes to create the
 * distributed synchronizer and distribute a mesh
 *
 * @section LICENSE
 *
 * Copyright (©) 2016-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "node_info_per_processor.hh"
#include "communicator.hh"
#include "node_group.hh"
#include "node_synchronizer.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
NodeInfoPerProc::NodeInfoPerProc(NodeSynchronizer & synchronizer,
                                 UInt message_cnt, UInt root)
    : MeshAccessor(synchronizer.getMesh()), synchronizer(synchronizer),
      comm(synchronizer.getCommunicator()), rank(comm.whoAmI()),
      nb_proc(comm.getNbProc()), root(root), mesh(synchronizer.getMesh()),
      spatial_dimension(synchronizer.getMesh().getSpatialDimension()),
      message_count(message_cnt) {}

/* -------------------------------------------------------------------------- */
void NodeInfoPerProc::synchronize() {
  synchronizeNodes();
  synchronizeTypes();
  synchronizeGroups();
  synchronizePeriodicity();
}

/* -------------------------------------------------------------------------- */
template <class CommunicationBuffer>
void NodeInfoPerProc::fillNodeGroupsFromBuffer(CommunicationBuffer & buffer) {
  AKANTU_DEBUG_IN();

  std::vector<std::vector<std::string>> node_to_group;

  buffer >> node_to_group;

  AKANTU_DEBUG_ASSERT(node_to_group.size() == mesh.getNbGlobalNodes(),
                      "Not the good amount of nodes where transmitted");

  const auto & global_nodes = mesh.getGlobalNodesIds();

  auto nbegin = global_nodes.begin();
  auto nit = global_nodes.begin();
  auto nend = global_nodes.end();
  for (; nit != nend; ++nit) {
    auto it = node_to_group[*nit].begin();
    auto end = node_to_group[*nit].end();

    for (; it != end; ++it) {
      mesh.getNodeGroup(*it).add(nit - nbegin, false);
    }
  }

  GroupManager::const_node_group_iterator ngi = mesh.node_group_begin();
  GroupManager::const_node_group_iterator nge = mesh.node_group_end();
  for (; ngi != nge; ++ngi) {
    NodeGroup & ng = *(ngi->second);
    ng.optimize();
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NodeInfoPerProc::fillNodesType() {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = mesh.getNbNodes();
  auto & nodes_flags = this->getNodesFlags();

  Array<UInt> nodes_set(nb_nodes);
  nodes_set.set(0);

  enum NodeSet {
    NORMAL_SET = 1,
    GHOST_SET = 2,
  };

  Array<bool> already_seen(nb_nodes, 1, false);

  for (auto gt : ghost_types) {
    UInt set = NORMAL_SET;
    if (gt == _ghost)
      set = GHOST_SET;

    already_seen.set(false);
    for (auto && type :
         mesh.elementTypes(_all_dimensions, gt, _ek_not_defined)) {
      const auto & connectivity = mesh.getConnectivity(type, gt);

      for (auto & conn :
           make_view(connectivity, connectivity.getNbComponent())) {
        for (UInt n = 0; n < conn.size(); ++n) {
          AKANTU_DEBUG_ASSERT(conn(n) < nb_nodes,
                              "Node " << conn(n)
                                      << " bigger than number of nodes "
                                      << nb_nodes);
          if (!already_seen(conn(n))) {
            nodes_set(conn(n)) += set;
            already_seen(conn(n)) = true;
          }
        }
      }
    }
  }

  nodes_flags.resize(nb_nodes);
  for (UInt i = 0; i < nb_nodes; ++i) {
    if (nodes_set(i) == NORMAL_SET)
      nodes_flags(i) = NodeFlag::_normal;
    else if (nodes_set(i) == GHOST_SET)
      nodes_flags(i) = NodeFlag::_pure_ghost;
    else if (nodes_set(i) == (GHOST_SET + NORMAL_SET))
      nodes_flags(i) = NodeFlag::_master;
    else {
      AKANTU_EXCEPTION("Gni ?");
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NodeInfoPerProc::fillCommunicationScheme(const Array<UInt> & master_info) {
  AKANTU_DEBUG_IN();

  Communications<UInt> & communications =
      this->synchronizer.getCommunications();

  { // send schemes
    auto it = master_info.begin_reinterpret(2, master_info.size() / 2);
    auto end = master_info.end_reinterpret(2, master_info.size() / 2);

    std::map<UInt, Array<UInt>> send_array_per_proc;

    for (; it != end; ++it) {
      const Vector<UInt> & send_info = *it;

      send_array_per_proc[send_info(0)].push_back(send_info(1));
    }

    for (auto & send_schemes : send_array_per_proc) {
      auto & scheme = communications.createSendScheme(send_schemes.first);

      auto & sends = send_schemes.second;
      std::sort(sends.begin(), sends.end());
      std::transform(sends.begin(), sends.end(), sends.begin(),
                     [this](UInt g) -> UInt { return mesh.getNodeLocalId(g); });
      scheme.copy(sends);
    }
  }

  { // receive schemes
    std::map<UInt, Array<UInt>> recv_array_per_proc;

    for (auto node : arange(mesh.getNbNodes())) {
      if (mesh.isSlaveNode(node)) {
        recv_array_per_proc[mesh.getNodePrank(node)].push_back(
            mesh.getNodeGlobalId(node));
      }
    }

    for (auto & recv_schemes : recv_array_per_proc) {
      auto & scheme = communications.createRecvScheme(recv_schemes.first);

      auto & recvs = recv_schemes.second;
      std::sort(recvs.begin(), recvs.end());
      std::transform(recvs.begin(), recvs.end(), recvs.begin(),
                     [this](UInt g) -> UInt { return mesh.getNodeLocalId(g); });

      scheme.copy(recvs);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NodeInfoPerProc::fillPeriodicPairs(const Array<UInt> & global_pairs,
                                        std::vector<UInt> & missing_nodes) {
  this->wipePeriodicInfo();
  auto & nodes_flags = this->getNodesFlags();

  auto checkIsLocal = [&](auto && global_node) {
    auto && node = mesh.getNodeLocalId(global_node);
    if (node == UInt(-1)) {
      auto & global_nodes = this->getNodesGlobalIds();
      node = global_nodes.size();

      global_nodes.push_back(global_node);
      nodes_flags.push_back(NodeFlag::_pure_ghost);
      missing_nodes.push_back(global_node);
      std::cout << "Missing node " << node << std::endl;
    }
    return node;
  };

  for (auto && pairs : make_view(global_pairs, 2)) {
    UInt slave = checkIsLocal(pairs(0));
    UInt master = checkIsLocal(pairs(1));

    this->addPeriodicSlave(slave, master);
  }

  this->markMeshPeriodic();
}

/* -------------------------------------------------------------------------- */
void NodeInfoPerProc::receiveMissingPeriodic(
    DynamicCommunicationBuffer & buffer) {
  auto & nodes = this->getNodes();
  Communications<UInt> & communications =
      this->synchronizer.getCommunications();

  std::size_t nb_nodes;
  buffer >> nb_nodes;

  for (auto _ [[gnu::unused]] : arange(nb_nodes)) {
    Vector<Real> pos(spatial_dimension);
    Int prank;
    buffer >> pos;
    buffer >> prank;

    UInt node = nodes.size();

    this->setNodePrank(node, prank);
    nodes.push_back(pos);

    auto & scheme = communications.createRecvScheme(prank);
    scheme.push_back(node);
  }

  while (buffer.getLeftToUnpack() != 0) {
    Int prank;
    UInt gnode;

    buffer >> gnode;
    buffer >> prank;
    auto node = mesh.getNodeLocalId(gnode);

    AKANTU_DEBUG_ASSERT(node != UInt(-1),
                        "I cannot send the node "
                            << gnode << " to proc " << prank
                            << " because it is note a local node");

    auto & scheme = communications.createSendScheme(prank);
    scheme.push_back(node);
  }
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
MasterNodeInfoPerProc::MasterNodeInfoPerProc(NodeSynchronizer & synchronizer,
                                             UInt message_cnt, UInt root)
    : NodeInfoPerProc(synchronizer, message_cnt, root),
      all_nodes(0, synchronizer.getMesh().getSpatialDimension()) {
  UInt nb_global_nodes = this->mesh.getNbGlobalNodes();
  this->comm.broadcast(nb_global_nodes, this->root);
}

/* -------------------------------------------------------------------------- */
void MasterNodeInfoPerProc::synchronizeNodes() {
  this->nodes_per_proc.resize(nb_proc);
  this->nb_nodes_per_proc.resize(nb_proc);

  Array<Real> local_nodes(0, spatial_dimension);
  Array<Real> & nodes = this->getNodes();

  all_nodes.copy(nodes);
  nodes_pranks.resize(nodes.size(), root);

  for (UInt p = 0; p < nb_proc; ++p) {
    UInt nb_nodes = 0;
    //      UInt * buffer;
    Array<Real> * nodes_to_send;

    Array<UInt> & nodespp = nodes_per_proc[p];
    if (p != root) {
      nodes_to_send = new Array<Real>(0, spatial_dimension);
      AKANTU_DEBUG_INFO("Receiving number of nodes from proc "
                        << p << " " << Tag::genTag(p, 0, Tag::_NB_NODES));
      comm.receive(nb_nodes, p, Tag::genTag(p, 0, Tag::_NB_NODES));
      nodespp.resize(nb_nodes);
      this->nb_nodes_per_proc(p) = nb_nodes;

      AKANTU_DEBUG_INFO("Receiving list of nodes from proc "
                        << p << " " << Tag::genTag(p, 0, Tag::_NODES));
      comm.receive(nodespp, p, Tag::genTag(p, 0, Tag::_NODES));
    } else {
      Array<UInt> & local_ids = this->getNodesGlobalIds();
      this->nb_nodes_per_proc(p) = local_ids.size();
      nodespp.copy(local_ids);
      nodes_to_send = &local_nodes;
    }

    Array<UInt>::const_scalar_iterator it = nodespp.begin();
    Array<UInt>::const_scalar_iterator end = nodespp.end();
    /// get the coordinates for the selected nodes
    for (; it != end; ++it) {
      Vector<Real> coord(nodes.storage() + spatial_dimension * *it,
                         spatial_dimension);
      nodes_to_send->push_back(coord);
    }

    if (p != root) { /// send them for distant processors
      AKANTU_DEBUG_INFO("Sending coordinates to proc "
                        << p << " "
                        << Tag::genTag(this->rank, 0, Tag::_COORDINATES));
      comm.send(*nodes_to_send, p,
                Tag::genTag(this->rank, 0, Tag::_COORDINATES));
      delete nodes_to_send;
    }
  }

  /// construct the local nodes coordinates
  nodes.copy(local_nodes);
}

/* -------------------------------------------------------------------------- */
void MasterNodeInfoPerProc::synchronizeTypes() {
  //           <global_id,     <proc, local_id> >
  std::multimap<UInt, std::pair<UInt, UInt>> nodes_to_proc;
  std::vector<Array<NodeFlag>> nodes_flags_per_proc(nb_proc);
  std::vector<Array<Int>> nodes_prank_per_proc(nb_proc);

  if (mesh.isPeriodic())
    all_periodic_flags.copy(this->getNodesFlags());

  // arrays containing pairs of (proc, node)
  std::vector<Array<UInt>> nodes_to_send_per_proc(nb_proc);
  for (UInt p = 0; p < nb_proc; ++p) {
    nodes_flags_per_proc[p].resize(nb_nodes_per_proc(p));
    nodes_prank_per_proc[p].resize(nb_nodes_per_proc(p));
  }

  this->fillNodesType();

  for (UInt p = 0; p < nb_proc; ++p) {
    auto & nodes_flags = nodes_flags_per_proc[p];
    if (p != root) {
      AKANTU_DEBUG_INFO(
          "Receiving first nodes types from proc "
          << p << " "
          << Tag::genTag(this->rank, this->message_count, Tag::_NODES_TYPE));
      comm.receive(nodes_flags, p, Tag::genTag(p, 0, Tag::_NODES_TYPE));
    } else {
      nodes_flags.copy(this->getNodesFlags());
    }

    // stack all processors claiming to be master for a node
    for (UInt local_node = 0; local_node < nb_nodes_per_proc(p); ++local_node) {
      if ((nodes_flags(local_node) & NodeFlag::_shared_mask) ==
          NodeFlag::_master) {
        UInt global_node = nodes_per_proc[p](local_node);
        nodes_to_proc.insert(
            std::make_pair(global_node, std::make_pair(p, local_node)));
      }
    }
  }

  for (UInt i = 0; i < mesh.getNbGlobalNodes(); ++i) {
    auto it_range = nodes_to_proc.equal_range(i);
    if (it_range.first == nodes_to_proc.end() || it_range.first->first != i)
      continue;

    // pick the first processor out of the multi-map as the actual master
    UInt master_proc = (it_range.first)->second.first;
    nodes_pranks[i] = master_proc;

    for (auto it_node = it_range.first; it_node != it_range.second; ++it_node) {
      UInt proc = it_node->second.first;
      UInt node = it_node->second.second;
      if (proc != master_proc) {
        // store the info on all the slaves for a given master
        nodes_flags_per_proc[proc](node) = NodeFlag::_slave;
        nodes_prank_per_proc[proc](node) = master_proc;
        nodes_to_send_per_proc[master_proc].push_back(proc);
        nodes_to_send_per_proc[master_proc].push_back(i);
      }
    }
  }

  std::vector<CommunicationRequest> requests_send_type;
  std::vector<CommunicationRequest> requests_send_master_info;
  for (UInt p = 0; p < nb_proc; ++p) {
    if (p != root) {
      AKANTU_DEBUG_INFO("Sending nodes types to proc "
                        << p << " "
                        << Tag::genTag(this->rank, 0, Tag::_NODES_TYPE));
      requests_send_type.push_back(
          comm.asyncSend(nodes_flags_per_proc[p], p,
                         Tag::genTag(this->rank, 0, Tag::_NODES_TYPE)));

      requests_send_type.push_back(
          comm.asyncSend(nodes_prank_per_proc[p], p,
                         Tag::genTag(this->rank, 2, Tag::_NODES_TYPE)));

      auto & nodes_to_send = nodes_to_send_per_proc[p];

      AKANTU_DEBUG_INFO("Sending nodes master info to proc "
                        << p << " "
                        << Tag::genTag(this->rank, 1, Tag::_NODES_TYPE));
      requests_send_master_info.push_back(comm.asyncSend(
          nodes_to_send, p, Tag::genTag(this->rank, 1, Tag::_NODES_TYPE)));
    } else {
      this->getNodesFlags().copy(nodes_flags_per_proc[p]);
      for (auto && data : enumerate(nodes_prank_per_proc[p])) {
        auto node = std::get<0>(data);
        if (mesh.isSlaveNode(node)) {
          this->setNodePrank(node, std::get<1>(data));
        }
      }

      this->fillCommunicationScheme(nodes_to_send_per_proc[root]);
    }
  }

  comm.waitAll(requests_send_type);
  comm.freeCommunicationRequest(requests_send_type);
  requests_send_type.clear();

  comm.waitAll(requests_send_master_info);
  comm.freeCommunicationRequest(requests_send_master_info);
}

/* -------------------------------------------------------------------------- */
void MasterNodeInfoPerProc::synchronizeGroups() {
  AKANTU_DEBUG_IN();

  UInt nb_total_nodes = mesh.getNbGlobalNodes();

  DynamicCommunicationBuffer buffer;

  using NodeToGroup = std::vector<std::vector<std::string>>;
  NodeToGroup node_to_group;
  node_to_group.resize(nb_total_nodes);

  GroupManager::const_node_group_iterator ngi = mesh.node_group_begin();
  GroupManager::const_node_group_iterator nge = mesh.node_group_end();
  for (; ngi != nge; ++ngi) {
    NodeGroup & ng = *(ngi->second);

    std::string name = ngi->first;

    NodeGroup::const_node_iterator nit = ng.begin();
    NodeGroup::const_node_iterator nend = ng.end();
    for (; nit != nend; ++nit) {
      node_to_group[*nit].push_back(name);
    }

    nit = ng.begin();
    if (nit != nend)
      ng.empty();
  }

  buffer << node_to_group;

  std::vector<CommunicationRequest> requests;
  for (UInt p = 0; p < nb_proc; ++p) {
    if (p == this->rank)
      continue;
    AKANTU_DEBUG_INFO("Sending node groups to proc "
                      << p << " "
                      << Tag::genTag(this->rank, p, Tag::_NODE_GROUP));
    requests.push_back(comm.asyncSend(
        buffer, p, Tag::genTag(this->rank, p, Tag::_NODE_GROUP)));
  }

  this->fillNodeGroupsFromBuffer(buffer);

  comm.waitAll(requests);
  comm.freeCommunicationRequest(requests);
  requests.clear();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MasterNodeInfoPerProc::synchronizePeriodicity() {
  bool is_periodic = mesh.isPeriodic();
  comm.broadcast(is_periodic, root);

  if (not is_periodic)
    return;

  std::vector<CommunicationRequest> requests;
  std::vector<Array<UInt>> periodic_info_to_send_per_proc;
  for (auto p : arange(nb_proc)) {
    periodic_info_to_send_per_proc.emplace_back(0, 2);
    auto && periodic_info = periodic_info_to_send_per_proc.back();

    for (UInt proc_local_node : arange(nb_nodes_per_proc(p))) {
      UInt global_node = nodes_per_proc[p](proc_local_node);
      if ((all_periodic_flags[global_node] & NodeFlag::_periodic_mask) ==
          NodeFlag::_periodic_slave) {
        periodic_info.push_back(
            Vector<UInt>{global_node, mesh.getPeriodicMaster(global_node)});
      }
    }

    if (p == root)
      continue;

    auto && tag = Tag::genTag(this->rank, p, Tag::_PERIODIC_SLAVES);
    AKANTU_DEBUG_INFO("Sending periodic info to proc " << p << " " << tag);
    requests.push_back(comm.asyncSend(periodic_info, p, tag));
  }

  CommunicationStatus status;
  std::vector<DynamicCommunicationBuffer> buffers(nb_proc);
  std::vector<std::vector<UInt>> proc_missings(nb_proc);
  auto nodes_it = all_nodes.begin(spatial_dimension);

  for (UInt p = 0; p < nb_proc; ++p) {
    auto & proc_missing = proc_missings[p];
    if (p != root) {
      auto && tag = Tag::genTag(p, 0, Tag::_PERIODIC_NODES);
      comm.probe<UInt>(p, tag, status);

      proc_missing.resize(status.size());
      comm.receive(proc_missing, p, tag);
    } else {
      fillPeriodicPairs(periodic_info_to_send_per_proc[root], proc_missing);
    }

    auto & buffer = buffers[p];
    buffer.reserve((spatial_dimension * sizeof(Real) + sizeof(Int)) *
                   proc_missing.size());
    buffer << proc_missing.size();
    for (auto && node : proc_missing) {
      buffer << *(nodes_it + node);
      buffer << nodes_pranks(node);
    }
  }

  for (UInt p = 0; p < nb_proc; ++p) {
    for (auto && node : proc_missings[p]) {
      auto & buffer = buffers[nodes_pranks(node)];
      buffer << node;
      buffer << p;
    }
  }

  for (UInt p = 0; p < nb_proc; ++p) {
    if (p != root) {
      auto && tag_send = Tag::genTag(p, 1, Tag::_PERIODIC_NODES);
      requests.push_back(comm.asyncSend(buffers[p], p, tag_send));
    } else {
      receiveMissingPeriodic(buffers[p]);
    }
  }

  comm.waitAll(requests);
  comm.freeCommunicationRequest(requests);
  requests.clear();
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
SlaveNodeInfoPerProc::SlaveNodeInfoPerProc(NodeSynchronizer & synchronizer,
                                           UInt message_cnt, UInt root)
    : NodeInfoPerProc(synchronizer, message_cnt, root) {
  UInt nb_global_nodes = 0;
  comm.broadcast(nb_global_nodes, root);
  this->setNbGlobalNodes(nb_global_nodes);
}

/* -------------------------------------------------------------------------- */
void SlaveNodeInfoPerProc::synchronizeNodes() {
  AKANTU_DEBUG_INFO("Sending list of nodes to proc "
                    << root << " " << Tag::genTag(this->rank, 0, Tag::_NB_NODES)
                    << " " << Tag::genTag(this->rank, 0, Tag::_NODES));
  Array<UInt> & local_ids = this->getNodesGlobalIds();
  Array<Real> & nodes = this->getNodes();

  UInt nb_nodes = local_ids.size();
  comm.send(nb_nodes, root, Tag::genTag(this->rank, 0, Tag::_NB_NODES));
  comm.send(local_ids, root, Tag::genTag(this->rank, 0, Tag::_NODES));

  /* --------<<<<-COORDINATES---------------------------------------------- */
  nodes.resize(nb_nodes);
  AKANTU_DEBUG_INFO("Receiving coordinates from proc "
                    << root << " " << Tag::genTag(root, 0, Tag::_COORDINATES));
  comm.receive(nodes, root, Tag::genTag(root, 0, Tag::_COORDINATES));
}

/* -------------------------------------------------------------------------- */
void SlaveNodeInfoPerProc::synchronizeTypes() {
  this->fillNodesType();

  auto & nodes_types = this->getNodesFlags();

  AKANTU_DEBUG_INFO("Sending first nodes types to proc "
                    << root << ""
                    << Tag::genTag(this->rank, 0, Tag::_NODES_TYPE));
  comm.send(nodes_types, root, Tag::genTag(this->rank, 0, Tag::_NODES_TYPE));

  AKANTU_DEBUG_INFO("Receiving nodes types from proc "
                    << root << " " << Tag::genTag(root, 0, Tag::_NODES_TYPE));
  comm.receive(nodes_types, root, Tag::genTag(root, 0, Tag::_NODES_TYPE));

  Array<Int> nodes_prank(nodes_types.size());
  comm.receive(nodes_prank, root, Tag::genTag(root, 2, Tag::_NODES_TYPE));
  for (auto && data : enumerate(nodes_prank)) {
    auto node = std::get<0>(data);
    if (mesh.isSlaveNode(node)) {
      this->setNodePrank(node, std::get<1>(data));
    }
  }

  AKANTU_DEBUG_INFO("Receiving nodes master info from proc "
                    << root << " " << Tag::genTag(root, 1, Tag::_NODES_TYPE));
  CommunicationStatus status;
  comm.probe<UInt>(root, Tag::genTag(root, 1, Tag::_NODES_TYPE), status);

  Array<UInt> nodes_master_info(status.size());
  if (nodes_master_info.size() > 0)
    comm.receive(nodes_master_info, root,
                 Tag::genTag(root, 1, Tag::_NODES_TYPE));

  this->fillCommunicationScheme(nodes_master_info);
}

/* -------------------------------------------------------------------------- */
void SlaveNodeInfoPerProc::synchronizeGroups() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Receiving node groups from proc "
                    << root << " "
                    << Tag::genTag(root, this->rank, Tag::_NODE_GROUP));

  CommunicationStatus status;
  comm.probe<char>(root, Tag::genTag(root, this->rank, Tag::_NODE_GROUP),
                   status);

  CommunicationBuffer buffer(status.size());
  comm.receive(buffer, root, Tag::genTag(root, this->rank, Tag::_NODE_GROUP));

  this->fillNodeGroupsFromBuffer(buffer);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SlaveNodeInfoPerProc::synchronizePeriodicity() {
  bool is_periodic;
  comm.broadcast(is_periodic, root);

  if (not is_periodic)
    return;

  CommunicationStatus status;
  auto && tag = Tag::genTag(root, this->rank, Tag::_PERIODIC_SLAVES);
  comm.probe<UInt>(root, tag, status);

  Array<UInt> periodic_info(status.size() / 2, 2);
  comm.receive(periodic_info, root, tag);

  std::vector<UInt> proc_missing;
  fillPeriodicPairs(periodic_info, proc_missing);

  auto && tag_missing_request =
      Tag::genTag(this->rank, 0, Tag::_PERIODIC_NODES);
  comm.send(proc_missing, root, tag_missing_request);

  DynamicCommunicationBuffer buffer;
  auto && tag_missing = Tag::genTag(this->rank, 1, Tag::_PERIODIC_NODES);
  comm.receive(buffer, root, tag_missing);

  receiveMissingPeriodic(buffer);
}

} // akantu
