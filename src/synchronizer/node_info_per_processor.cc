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
  Array<NodeType> & nodes_type = this->getNodesType();

  Array<UInt> nodes_set(nb_nodes);
  nodes_set.set(0);

  enum NodeSet {
    NORMAL_SET = 1,
    GHOST_SET = 2,
  };

  Array<bool> already_seen(nb_nodes, 1, false);

  for (UInt g = _not_ghost; g <= _ghost; ++g) {
    auto gt = (GhostType)g;
    UInt set = NORMAL_SET;
    if (gt == _ghost)
      set = GHOST_SET;

    already_seen.set(false);
    Mesh::type_iterator it =
        mesh.firstType(_all_dimensions, gt, _ek_not_defined);
    Mesh::type_iterator end =
        mesh.lastType(_all_dimensions, gt, _ek_not_defined);
    for (; it != end; ++it) {
      ElementType type = *it;

      UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
      UInt nb_element = mesh.getNbElement(type, gt);
      Array<UInt>::const_vector_iterator conn_it =
          mesh.getConnectivity(type, gt).begin(nb_nodes_per_element);

      for (UInt e = 0; e < nb_element; ++e, ++conn_it) {
        const Vector<UInt> & conn = *conn_it;
        for (UInt n = 0; n < nb_nodes_per_element; ++n) {
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

  for (UInt i = 0; i < nb_nodes; ++i) {
    if (nodes_set(i) == NORMAL_SET)
      nodes_type(i) = _nt_normal;
    else if (nodes_set(i) == GHOST_SET)
      nodes_type(i) = _nt_pure_gost;
    else if (nodes_set(i) == (GHOST_SET + NORMAL_SET))
      nodes_type(i) = _nt_master;
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
    const Array<NodeType> & nodes_type = this->getNodesType();

    std::map<UInt, Array<UInt>> recv_array_per_proc;

    UInt node = 0;
    for (auto & node_type : nodes_type) {
      if (Int(node_type) >= 0) {
        recv_array_per_proc[node_type].push_back(mesh.getNodeGlobalId(node));
      }
      ++node;
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
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
MasterNodeInfoPerProc::MasterNodeInfoPerProc(NodeSynchronizer & synchronizer,
                                             UInt message_cnt, UInt root)
    : NodeInfoPerProc(synchronizer, message_cnt, root) {
  UInt nb_global_nodes = this->mesh.getNbGlobalNodes();
  this->comm.broadcast(nb_global_nodes, this->root);
}

/* -------------------------------------------------------------------------- */
void MasterNodeInfoPerProc::synchronizeNodes() {
  this->nodes_per_proc.resize(nb_proc);
  this->nb_nodes_per_proc.resize(nb_proc);

  Array<Real> local_nodes(0, spatial_dimension);
  Array<Real> & nodes = this->getNodes();

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
  std::vector<Array<NodeType>> nodes_type_per_proc(nb_proc);

  // arrays containing pairs of (proc, node)
  std::vector<Array<UInt>> nodes_to_send_per_proc(nb_proc);
  for (UInt p = 0; p < nb_proc; ++p) {
    nodes_type_per_proc[p].resize(nb_nodes_per_proc(p));
  }

  this->fillNodesType();

  for (UInt p = 0; p < nb_proc; ++p) {
    auto & nodes_types = nodes_type_per_proc[p];
    if (p != root) {
      AKANTU_DEBUG_INFO(
          "Receiving first nodes types from proc "
          << p << " "
          << Tag::genTag(this->rank, this->message_count, Tag::_NODES_TYPE));
      comm.receive(nodes_types, p, Tag::genTag(p, 0, Tag::_NODES_TYPE));
    } else {
      nodes_types.copy(this->getNodesType());
    }

    // stack all processors claiming to be master for a node
    for (UInt local_node = 0; local_node < nb_nodes_per_proc(p); ++local_node) {
      if (nodes_types(local_node) == _nt_master) {
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
    for (auto it_node = it_range.first; it_node != it_range.second; ++it_node) {
      UInt proc = it_node->second.first;
      UInt node = it_node->second.second;
      if (proc != master_proc) {
        // store the info on all the slaves for a given master
        nodes_type_per_proc[proc](node) = NodeType(master_proc);
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
          comm.asyncSend(nodes_type_per_proc[p], p,
                         Tag::genTag(this->rank, 0, Tag::_NODES_TYPE)));

      auto & nodes_to_send = nodes_to_send_per_proc[p];

      /// push back an element to avoid a send of size 0
      nodes_to_send.push_back(-1);

      AKANTU_DEBUG_INFO("Sending nodes master info to proc "
                        << p << " "
                        << Tag::genTag(this->rank, 1, Tag::_NODES_TYPE));
      requests_send_master_info.push_back(comm.asyncSend(
          nodes_to_send, p, Tag::genTag(this->rank, 1, Tag::_NODES_TYPE)));
    } else {
      this->getNodesType().copy(nodes_type_per_proc[p]);

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

  Array<NodeType> & nodes_types = this->getNodesType();

  AKANTU_DEBUG_INFO("Sending first nodes types to proc "
                    << root << ""
                    << Tag::genTag(this->rank, 0, Tag::_NODES_TYPE));
  comm.send(nodes_types, root, Tag::genTag(this->rank, 0, Tag::_NODES_TYPE));

  AKANTU_DEBUG_INFO("Receiving nodes types from proc "
                    << root << " " << Tag::genTag(root, 0, Tag::_NODES_TYPE));
  comm.receive(nodes_types, root, Tag::genTag(root, 0, Tag::_NODES_TYPE));

  AKANTU_DEBUG_INFO("Receiving nodes master info from proc "
                    << root << " " << Tag::genTag(root, 1, Tag::_NODES_TYPE));
  CommunicationStatus status;
  comm.probe<UInt>(root, Tag::genTag(root, 1, Tag::_NODES_TYPE), status);

  Array<UInt> nodes_master_info(status.size());

  comm.receive(nodes_master_info, root, Tag::genTag(root, 1, Tag::_NODES_TYPE));

  nodes_master_info.resize(nodes_master_info.size() - 1);

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

} // akantu
