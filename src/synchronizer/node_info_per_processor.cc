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
  synchronizeTags();
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

  for (auto && data : enumerate(global_nodes)) {
    for (const auto & node : node_to_group[std::get<1>(data)]) {
      mesh.getNodeGroup(node).add(std::get<0>(data), false);
    }
  }

  for (auto && ng_data : mesh.iterateNodeGroups()) {
    ng_data.optimize();
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
    if (gt == _ghost) {
      set = GHOST_SET;
    }

    already_seen.set(false);
    for (auto && type :
         mesh.elementTypes(_all_dimensions, gt, _ek_not_defined)) {
      const auto & connectivity = mesh.getConnectivity(type, gt);

      for (const auto & conn :
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
    if (nodes_set(i) == NORMAL_SET) {
      nodes_flags(i) = NodeFlag::_normal;
    } else if (nodes_set(i) == GHOST_SET) {
      nodes_flags(i) = NodeFlag::_pure_ghost;
    } else if (nodes_set(i) == (GHOST_SET + NORMAL_SET)) {
      nodes_flags(i) = NodeFlag::_master;
    } else {
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
    std::map<UInt, Array<UInt>> send_array_per_proc;

    for (const auto & send_info : make_view(master_info, 2)) {
      send_array_per_proc[send_info(0)].push_back(send_info(1));
    }

    for (auto & send_schemes : send_array_per_proc) {
      auto & scheme = communications.createSendScheme(send_schemes.first);

      auto & sends = send_schemes.second;
      std::sort(sends.begin(), sends.end());
      std::transform(sends.begin(), sends.end(), sends.begin(),
                     [this](UInt g) -> UInt { return mesh.getNodeLocalId(g); });
      scheme.copy(sends);
      AKANTU_DEBUG_INFO("Proc " << rank << " has " << sends.size()
                                << " nodes to send to  to proc "
                                << send_schemes.first);
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
      AKANTU_DEBUG_INFO("Proc " << rank << " will receive " << recvs.size()
                                << " nodes from proc " << recv_schemes.first);
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
void NodeInfoPerProc::fillNodalData(DynamicCommunicationBuffer & buffer,
                                    const std::string & tag_name) {

#define AKANTU_DISTRIBUTED_SYNHRONIZER_TAG_DATA(r, _, elem)                    \
  case MeshDataTypeCode::BOOST_PP_TUPLE_ELEM(2, 0, elem): {                    \
    auto & nodal_data =                                                        \
        mesh.getNodalData<BOOST_PP_TUPLE_ELEM(2, 1, elem)>(tag_name);          \
    nodal_data.resize(mesh.getNbNodes());                                      \
    for (auto && data : make_view(nodal_data)) {                               \
      buffer >> data;                                                          \
    }                                                                          \
    break;                                                                     \
  }

  MeshDataTypeCode data_type_code =
      mesh.getTypeCode(tag_name, MeshDataType::_nodal);
  switch (data_type_code) {
    BOOST_PP_SEQ_FOR_EACH(AKANTU_DISTRIBUTED_SYNHRONIZER_TAG_DATA, ,
                          AKANTU_MESH_DATA_TYPES)
  default:
    AKANTU_ERROR("Could not obtain the type of tag" << tag_name << "!");
    break;
  }
#undef AKANTU_DISTRIBUTED_SYNHRONIZER_TAG_DATA
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
  nodes_pranks.resize(nodes.size(), UInt(-1));

  for (UInt p = 0; p < nb_proc; ++p) {
    UInt nb_nodes = 0;
    //      UInt * buffer;
    Array<Real> * nodes_to_send{nullptr};

    Array<UInt> & nodespp = nodes_per_proc[p];
    if (p != root) {
      nodes_to_send = new Array<Real>(0, spatial_dimension);
      AKANTU_DEBUG_INFO("Receiving number of nodes from proc "
                        << p << " " << Tag::genTag(p, 0, Tag::_nb_nodes));
      comm.receive(nb_nodes, p, Tag::genTag(p, 0, Tag::_nb_nodes));
      nodespp.resize(nb_nodes);
      this->nb_nodes_per_proc(p) = nb_nodes;

      AKANTU_DEBUG_INFO("Receiving list of nodes from proc "
                        << p << " " << Tag::genTag(p, 0, Tag::_nodes));
      comm.receive(nodespp, p, Tag::genTag(p, 0, Tag::_nodes));
    } else {
      Array<UInt> & local_ids = this->getNodesGlobalIds();
      this->nb_nodes_per_proc(p) = local_ids.size();
      nodespp.copy(local_ids);
      nodes_to_send = &local_nodes;
    }

    /// get the coordinates for the selected nodes
    for (const auto & node : nodespp) {
      Vector<Real> coord(nodes.storage() + spatial_dimension * node,
                         spatial_dimension);
      nodes_to_send->push_back(coord);
    }

    if (p != root) { /// send them for distant processors
      AKANTU_DEBUG_INFO("Sending coordinates to proc "
                        << p << " "
                        << Tag::genTag(this->rank, 0, Tag::_coordinates));
      comm.send(*nodes_to_send, p,
                Tag::genTag(this->rank, 0, Tag::_coordinates));

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

  if (mesh.isPeriodic()) {
    all_periodic_flags.copy(this->getNodesFlags());
  }

  // arrays containing pairs of (proc, node)
  std::vector<Array<UInt>> nodes_to_send_per_proc(nb_proc);
  for (UInt p = 0; p < nb_proc; ++p) {
    nodes_flags_per_proc[p].resize(nb_nodes_per_proc(p), NodeFlag(0xFF));
    nodes_prank_per_proc[p].resize(nb_nodes_per_proc(p), -1);
  }

  this->fillNodesType();

  auto is_master = [](auto && flag) {
    return (flag & NodeFlag::_shared_mask) == NodeFlag::_master;
  };

  auto is_local = [](auto && flag) {
    return (flag & NodeFlag::_shared_mask) == NodeFlag::_normal;
  };

  for (auto p : arange(nb_proc)) {
    auto & nodes_flags = nodes_flags_per_proc[p];

    if (p != root) {
      AKANTU_DEBUG_INFO(
          "Receiving first nodes types from proc "
          << p << " "
          << Tag::genTag(this->rank, this->message_count, Tag::_nodes_type));
      comm.receive(nodes_flags, p, Tag::genTag(p, 0, Tag::_nodes_type));
    } else {
      nodes_flags.copy(this->getNodesFlags());
    }

    // stack all processors claiming to be master for a node
    for (auto local_node : arange(nb_nodes_per_proc(p))) {
      auto global_node = nodes_per_proc[p](local_node);

      if (is_master(nodes_flags(local_node))) {
        nodes_to_proc.insert(
            std::make_pair(global_node, std::make_pair(p, local_node)));

      } else if (is_local(nodes_flags(local_node))) {
        nodes_pranks[global_node] = p;
      }
    }
  }

  for (auto i : arange(mesh.getNbGlobalNodes())) {
    auto it_range = nodes_to_proc.equal_range(i);
    if (it_range.first == nodes_to_proc.end() || it_range.first->first != i) {
      continue;
    }

    // pick the first processor out of the multi-map as the actual master
    UInt master_proc = (it_range.first)->second.first;
    nodes_pranks[i] = master_proc;

    for (auto && data : range(it_range.first, it_range.second)) {
      auto proc = data.second.first;
      auto node = data.second.second;

      if (proc != master_proc) {
        // store the info on all the slaves for a given master
        nodes_flags_per_proc[proc](node) = NodeFlag::_slave;
        nodes_to_send_per_proc[master_proc].push_back(proc);
        nodes_to_send_per_proc[master_proc].push_back(i);
      }
    }
  }

  /// Fills the nodes prank per proc
  for (auto && data : zip(arange(nb_proc), nodes_per_proc, nodes_prank_per_proc,
                          nodes_flags_per_proc)) {
    for (auto && node_data :
         zip(std::get<1>(data), std::get<2>(data), std::get<3>(data))) {
      if (std::get<2>(node_data) == NodeFlag::_normal) {
        std::get<1>(node_data) = std::get<0>(data);
      } else {
        std::get<1>(node_data) = nodes_pranks(std::get<0>(node_data));
      }
    }
  }

  std::vector<CommunicationRequest> requests_send_type;
  std::vector<CommunicationRequest> requests_send_master_info;
  for (UInt p = 0; p < nb_proc; ++p) {
    if (p != root) {
      auto tag0 = Tag::genTag(this->rank, 0, Tag::_nodes_type);
      AKANTU_DEBUG_INFO("Sending nodes types to proc " << p << " " << tag0);
      requests_send_type.push_back(
          comm.asyncSend(nodes_flags_per_proc[p], p, tag0));

      auto tag2 = Tag::genTag(this->rank, 2, Tag::_nodes_type);
      AKANTU_DEBUG_INFO("Sending nodes pranks to proc " << p << " " << tag2);
      requests_send_type.push_back(
          comm.asyncSend(nodes_prank_per_proc[p], p, tag2));

      auto & nodes_to_send = nodes_to_send_per_proc[p];

      auto tag1 = Tag::genTag(this->rank, 1, Tag::_nodes_type);
      AKANTU_DEBUG_INFO("Sending nodes master info to proc " << p << " "
                                                             << tag1);
      requests_send_master_info.push_back(
          comm.asyncSend(nodes_to_send, p, tag1));
    } else {
      this->getNodesFlags().copy(nodes_flags_per_proc[p]);
      for (auto && data : enumerate(nodes_prank_per_proc[p])) {
        auto node = std::get<0>(data);
        if (not(mesh.isMasterNode(node) or mesh.isLocalNode(node))) {
          this->setNodePrank(node, std::get<1>(data));
        }
      }

      this->fillCommunicationScheme(nodes_to_send_per_proc[root]);
    }
  }

  Communicator::waitAll(requests_send_type);
  Communicator::freeCommunicationRequest(requests_send_type);

  Communicator::waitAll(requests_send_master_info);
  Communicator::freeCommunicationRequest(requests_send_master_info);
}

/* -------------------------------------------------------------------------- */
void MasterNodeInfoPerProc::synchronizeGroups() {
  AKANTU_DEBUG_IN();

  UInt nb_total_nodes = mesh.getNbGlobalNodes();

  DynamicCommunicationBuffer buffer;

  using NodeToGroup = std::vector<std::vector<std::string>>;
  NodeToGroup node_to_group;
  node_to_group.resize(nb_total_nodes);

  for (auto & ng : mesh.iterateNodeGroups()) {
    std::string name = ng.getName();

    for (auto && node : ng.getNodes()) {
      node_to_group[node].push_back(name);
    }

    ng.clear();
  }

  buffer << node_to_group;

  std::vector<CommunicationRequest> requests;
  for (UInt p = 0; p < nb_proc; ++p) {
    if (p == this->rank) {
      continue;
    }
    AKANTU_DEBUG_INFO("Sending node groups to proc "
                      << p << " "
                      << Tag::genTag(this->rank, p, Tag::_node_group));
    requests.push_back(comm.asyncSend(
        buffer, p, Tag::genTag(this->rank, p, Tag::_node_group)));
  }

  this->fillNodeGroupsFromBuffer(buffer);

  Communicator::waitAll(requests);
  Communicator::freeCommunicationRequest(requests);
  requests.clear();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MasterNodeInfoPerProc::synchronizePeriodicity() {
  bool is_periodic = mesh.isPeriodic();
  comm.broadcast(is_periodic, root);

  if (not is_periodic) {
    return;
  }

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

    if (p == root) {
      continue;
    }

    auto && tag = Tag::genTag(this->rank, p, Tag::_periodic_slaves);
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
      auto && tag = Tag::genTag(p, 0, Tag::_periodic_nodes);
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
      auto && tag_send = Tag::genTag(p, 1, Tag::_periodic_nodes);
      requests.push_back(comm.asyncSend(buffers[p], p, tag_send));
    } else {
      receiveMissingPeriodic(buffers[p]);
    }
  }

  Communicator::waitAll(requests);
  Communicator::freeCommunicationRequest(requests);
}

/* -------------------------------------------------------------------------- */
void MasterNodeInfoPerProc::fillTagBuffers(
    std::vector<DynamicCommunicationBuffer> & buffers,
    const std::string & tag_name) {
#define AKANTU_DISTRIBUTED_SYNHRONIZER_TAG_DATA(r, _, elem)                    \
  case MeshDataTypeCode::BOOST_PP_TUPLE_ELEM(2, 0, elem): {                    \
    auto & nodal_data =                                                        \
        mesh.getNodalData<BOOST_PP_TUPLE_ELEM(2, 1, elem)>(tag_name);          \
    for (auto && data : enumerate(nodes_per_proc)) {                           \
      auto proc = std::get<0>(data);                                           \
      auto & nodes = std::get<1>(data);                                        \
      auto & buffer = buffers[proc];                                           \
      for (auto & node : nodes) {                                              \
        for (auto i : arange(nodal_data.getNbComponent())) {                   \
          buffer << nodal_data(node, i);                                       \
        }                                                                      \
      }                                                                        \
    }                                                                          \
    break;                                                                     \
  }

  MeshDataTypeCode data_type_code =
      mesh.getTypeCode(tag_name, MeshDataType::_nodal);
  switch (data_type_code) {
    BOOST_PP_SEQ_FOR_EACH(AKANTU_DISTRIBUTED_SYNHRONIZER_TAG_DATA, ,
                          AKANTU_MESH_DATA_TYPES)
  default:
    AKANTU_ERROR("Could not obtain the type of tag" << tag_name << "!");
    break;
  }
#undef AKANTU_DISTRIBUTED_SYNHRONIZER_TAG_DATA
} // namespace akantu

/* -------------------------------------------------------------------------- */
void MasterNodeInfoPerProc::synchronizeTags() {
  /// tag info
  auto tag_names = mesh.getTagNames();

  DynamicCommunicationBuffer tags_buffer;
  for (auto && tag_name : tag_names) {
    tags_buffer << tag_name;
    tags_buffer << mesh.getTypeCode(tag_name, MeshDataType::_nodal);
    tags_buffer << mesh.getNbComponent(tag_name);
  }

  AKANTU_DEBUG_INFO(
      "Broadcasting the information about the nodes mesh data tags: ("
      << tags_buffer.size() << ").");
  comm.broadcast(tags_buffer, root);

  for (auto && tag_data : enumerate(tag_names)) {
    auto tag_count = std::get<0>(tag_data);
    auto & tag_name = std::get<1>(tag_data);
    std::vector<DynamicCommunicationBuffer> buffers;
    std::vector<CommunicationRequest> requests;
    buffers.resize(nb_proc);
    fillTagBuffers(buffers, tag_name);
    for (auto && data : enumerate(buffers)) {
      auto && proc = std::get<0>(data);
      auto & buffer = std::get<1>(data);
      if (proc == root) {
        fillNodalData(buffer, tag_name);
      } else {
        auto && tag = Tag::genTag(this->rank, tag_count, Tag::_mesh_data);
        requests.push_back(comm.asyncSend(buffer, proc, tag));
      }
    }

    Communicator::waitAll(requests);
    Communicator::freeCommunicationRequest(requests);
  }
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
                    << root << " " << Tag::genTag(this->rank, 0, Tag::_nb_nodes)
                    << " " << Tag::genTag(this->rank, 0, Tag::_nodes));
  Array<UInt> & local_ids = this->getNodesGlobalIds();
  Array<Real> & nodes = this->getNodes();

  UInt nb_nodes = local_ids.size();
  comm.send(nb_nodes, root, Tag::genTag(this->rank, 0, Tag::_nb_nodes));
  comm.send(local_ids, root, Tag::genTag(this->rank, 0, Tag::_nodes));

  /* --------<<<<-COORDINATES---------------------------------------------- */
  nodes.resize(nb_nodes);
  AKANTU_DEBUG_INFO("Receiving coordinates from proc "
                    << root << " " << Tag::genTag(root, 0, Tag::_coordinates));
  comm.receive(nodes, root, Tag::genTag(root, 0, Tag::_coordinates));
}

/* -------------------------------------------------------------------------- */
void SlaveNodeInfoPerProc::synchronizeTypes() {
  this->fillNodesType();

  auto & nodes_flags = this->getNodesFlags();

  AKANTU_DEBUG_INFO("Sending first nodes types to proc "
                    << root << ""
                    << Tag::genTag(this->rank, 0, Tag::_nodes_type));
  comm.send(nodes_flags, root, Tag::genTag(this->rank, 0, Tag::_nodes_type));

  AKANTU_DEBUG_INFO("Receiving nodes types from proc "
                    << root << " " << Tag::genTag(root, 0, Tag::_nodes_type));
  comm.receive(nodes_flags, root, Tag::genTag(root, 0, Tag::_nodes_type));

  Array<Int> nodes_prank(nodes_flags.size());

  AKANTU_DEBUG_INFO("Receiving nodes pranks from proc "
                    << root << " " << Tag::genTag(root, 2, Tag::_nodes_type));
  comm.receive(nodes_prank, root, Tag::genTag(root, 2, Tag::_nodes_type));
  for (auto && data : enumerate(nodes_prank)) {
    auto node = std::get<0>(data);
    if (not(mesh.isMasterNode(node) or mesh.isLocalNode(node))) {
      this->setNodePrank(node, std::get<1>(data));
    }
  }

  AKANTU_DEBUG_INFO("Receiving nodes master info from proc "
                    << root << " " << Tag::genTag(root, 1, Tag::_nodes_type));
  CommunicationStatus status;
  comm.probe<UInt>(root, Tag::genTag(root, 1, Tag::_nodes_type), status);

  Array<UInt> nodes_master_info(status.size());
  comm.receive(nodes_master_info, root, Tag::genTag(root, 1, Tag::_nodes_type));

  this->fillCommunicationScheme(nodes_master_info);
}

/* -------------------------------------------------------------------------- */
void SlaveNodeInfoPerProc::synchronizeGroups() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Receiving node groups from proc "
                    << root << " "
                    << Tag::genTag(root, this->rank, Tag::_node_group));

  DynamicCommunicationBuffer buffer;
  comm.receive(buffer, root, Tag::genTag(root, this->rank, Tag::_node_group));

  this->fillNodeGroupsFromBuffer(buffer);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SlaveNodeInfoPerProc::synchronizePeriodicity() {
  bool is_periodic;
  comm.broadcast(is_periodic, root);

  if (not is_periodic) {
    return;
  }

  CommunicationStatus status;
  auto && tag = Tag::genTag(root, this->rank, Tag::_periodic_slaves);
  comm.probe<UInt>(root, tag, status);

  Array<UInt> periodic_info(status.size() / 2, 2);
  comm.receive(periodic_info, root, tag);

  std::vector<UInt> proc_missing;
  fillPeriodicPairs(periodic_info, proc_missing);

  auto && tag_missing_request =
      Tag::genTag(this->rank, 0, Tag::_periodic_nodes);
  comm.send(proc_missing, root, tag_missing_request);

  DynamicCommunicationBuffer buffer;
  auto && tag_missing = Tag::genTag(this->rank, 1, Tag::_periodic_nodes);
  comm.receive(buffer, root, tag_missing);

  receiveMissingPeriodic(buffer);
}

/* -------------------------------------------------------------------------- */
void SlaveNodeInfoPerProc::synchronizeTags() {
  DynamicCommunicationBuffer tags_buffer;
  comm.broadcast(tags_buffer, root);

  std::vector<std::string> tag_names;
  while (tags_buffer.getLeftToUnpack() > 0) {
    std::string name;
    MeshDataTypeCode code;
    UInt nb_components;
    tags_buffer >> name;
    tags_buffer >> code;
    tags_buffer >> nb_components;

    mesh.registerNodalData(name, nb_components, code);
    tag_names.push_back(name);
  }

  for (auto && tag_data : enumerate(tag_names)) {
    auto tag_count = std::get<0>(tag_data);
    auto & tag_name = std::get<1>(tag_data);
    DynamicCommunicationBuffer buffer;
    auto && tag = Tag::genTag(this->root, tag_count, Tag::_mesh_data);
    comm.receive(buffer, this->root, tag);

    fillNodalData(buffer, tag_name);
  }
}

} // namespace akantu
