/**
 * @file   mesh_periodic.cc
 *
 * @author Nicolas Richart
 *
 * @date creation  Sat Feb 10 2018
 *
 * @brief Implementation of the perdiodicity capabilities in the mesh
 *
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "communication_tag.hh"
#include "communicator.hh"
#include "element_group.hh"
#include "mesh.hh"
#include "periodic_node_synchronizer.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
void Mesh::makePeriodic(const SpatialDirection & direction) {
  Array<UInt> list_1;
  Array<UInt> list_2;
  Real tolerance = 1e-10;

  auto lower_bound = this->getLowerBounds();
  auto upper_bound = this->getUpperBounds();
  auto length = upper_bound(direction) - lower_bound(direction);

  const auto & positions = *nodes;

  for (auto && data : enumerate(make_view(positions, spatial_dimension))) {
    UInt node = std::get<0>(data);
    const auto & pos = std::get<1>(data);

    if (std::abs((pos(direction) - lower_bound(direction)) / length) <
        tolerance) {
      list_1.push_back(node);
    }

    if (std::abs((pos(direction) - upper_bound(direction)) / length) <
        tolerance) {
      list_2.push_back(node);
    }
  }

  this->makePeriodic(direction, list_1, list_2);
}

/* -------------------------------------------------------------------------- */
void Mesh::makePeriodic(const SpatialDirection & direction, const ID & list_1,
                        const ID & list_2) {
  const auto & list_nodes_1 =
      mesh.getElementGroup(list_1).getNodeGroup().getNodes();
  const auto & list_nodes_2 =
      mesh.getElementGroup(list_2).getNodeGroup().getNodes();

  this->makePeriodic(direction, list_nodes_1, list_nodes_2);
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

namespace {
  struct NodeInfo {
    NodeInfo() = default;
    NodeInfo(UInt spatial_dimension) : position(spatial_dimension) {}
    NodeInfo(UInt node, const Vector<Real> & position,
             const SpatialDirection & direction)
        : node(node), position(position) {
      this->direction_position = position(direction);
      this->position(direction) = 0.;
    }

    NodeInfo(const NodeInfo & other) = default;
    NodeInfo(NodeInfo && other) noexcept = default;
    NodeInfo & operator=(const NodeInfo & other) = default;
    NodeInfo & operator=(NodeInfo && other)  = default;

    UInt node{0};
    Vector<Real> position;
    Real direction_position{0.};
  };
} // namespace

/* -------------------------------------------------------------------------- */
// left is for lower values on direction and right for highest values
void Mesh::makePeriodic(const SpatialDirection & direction,
                        const Array<UInt> & list_left,
                        const Array<UInt> & list_right) {
  Real tolerance = 1e-10;

  const auto & positions = *nodes;
  auto lower_bound = this->getLowerBounds();
  auto upper_bound = this->getUpperBounds();
  auto length = upper_bound(direction) - lower_bound(direction);

  lower_bound(direction) = 0;
  upper_bound(direction) = 0;

  auto prank = communicator->whoAmI();

  std::vector<NodeInfo> nodes_left(list_left.size());
  std::vector<NodeInfo> nodes_right(list_right.size());

  BBox bbox(spatial_dimension);
  auto to_position = [&](UInt node) {
    Vector<Real> pos(spatial_dimension);
    for (UInt s : arange(spatial_dimension)) {
      pos(s) = positions(node, s);
    }
    auto && info = NodeInfo(node, pos, direction);
    bbox += info.position;
    return std::move(info);
  };

  std::transform(list_left.begin(), list_left.end(), nodes_left.begin(),
                 to_position);
  BBox bbox_left = bbox;

  bbox.reset();
  std::transform(list_right.begin(), list_right.end(), nodes_right.begin(),
                 to_position);
  BBox bbox_right = bbox;

  std::vector<UInt> new_nodes;
  if (is_distributed) {
    NewNodesEvent event(AKANTU_CURRENT_FUNCTION);

    /* ---------------------------------------------------------------------- */
    // function to send nodes in bboxes intersections
    auto extract_and_send_nodes = [&](const auto & bbox, const auto & node_list,
                                      auto & buffers, auto proc, auto cnt) {
      // buffers.resize(buffers.size() + 1);
      buffers.push_back(std::make_unique<DynamicCommunicationBuffer>());
      auto & buffer = *buffers.back();

      // std::cout << "Sending to " << proc << std::endl;
      for (auto & info : node_list) {
        if (bbox.contains(info.position) and isLocalOrMasterNode(info.node)) {
          Vector<Real> pos = info.position;
          pos(direction) = info.direction_position;

          NodeFlag flag = (*nodes_flags)(info.node) & NodeFlag::_periodic_mask;
          UInt gnode = getNodeGlobalId(info.node);
          buffer << gnode;
          buffer << pos;
          buffer << flag;

          // std::cout << " - node " << getNodeGlobalId(info.node);

          // if is slave sends master info
          if (flag == NodeFlag::_periodic_slave) {
            UInt master = getNodeGlobalId(periodic_slave_master[info.node]);
            // std::cout << " slave of " << master << std::endl;
            buffer << master;
          }

          // if is master sends list of slaves
          if (flag == NodeFlag::_periodic_master) {
            UInt nb_slaves = periodic_master_slave.count(info.node);
            buffer << nb_slaves;

            // std::cout << " master of " << nb_slaves << " nodes : [";
            auto slaves = periodic_master_slave.equal_range(info.node);
            for (auto it = slaves.first; it != slaves.second; ++it) {
              UInt gslave = getNodeGlobalId(it->second);
              // std::cout << (it == slaves.first ? "" : ", ") << gslave;
              buffer << gslave;
            }
            // std::cout << "]";
          }
          // std::cout << std::endl;
        }
      }

      auto tag = Tag::genTag(prank, 10 * direction + cnt, Tag::_periodic_nodes);
      // std::cout << "SBuffer size " << buffer.size() << " " << tag <<
      // std::endl;
      return communicator->asyncSend(buffer, proc, tag);
    };

    /* ---------------------------------------------------------------------- */
    // function to receive nodes in bboxes intersections
    auto recv_and_extract_nodes = [&](auto & node_list, const auto proc,
                                      auto cnt) {
      DynamicCommunicationBuffer buffer;
      auto tag = Tag::genTag(proc, 10 * direction + cnt, Tag::_periodic_nodes);
      communicator->receive(buffer, proc, tag);
      // std::cout << "RBuffer size " << buffer.size() << " " << tag <<
      // std::endl; std::cout << "Receiving from " << proc << std::endl;

      while (not buffer.empty()) {
        Vector<Real> pos(spatial_dimension);
        UInt global_node;
        NodeFlag flag;
        buffer >> global_node;
        buffer >> pos;
        buffer >> flag;

        // std::cout << " - node " << global_node;
        auto local_node = getNodeLocalId(global_node);

        // get the master info of is slave
        if (flag == NodeFlag::_periodic_slave) {
          UInt master_node;
          buffer >> master_node;
          // std::cout << " slave of " << master_node << std::endl;
          // auto local_master_node = getNodeLocalId(master_node);
          // AKANTU_DEBUG_ASSERT(local_master_node != UInt(-1),
          //"Should I know the master node " << master_node);
        }

        // get the list of slaves if is master
        if ((flag & NodeFlag::_periodic_mask) == NodeFlag::_periodic_master) {
          UInt nb_slaves;
          buffer >> nb_slaves;
          // std::cout << " master of " << nb_slaves << " nodes : [";
          for (auto ns [[gnu::unused]] : arange(nb_slaves)) {
            UInt gslave_node;
            buffer >> gslave_node;
            // std::cout << (ns == 0 ? "" : ", ") << gslave_node;
            // auto lslave_node = getNodeLocalId(gslave_node);
            // AKANTU_DEBUG_ASSERT(lslave_node != UInt(-1),
            //                    "Should I know the slave node " <<
            //                    gslave_node);
          }
          // std::cout << "]";
        }
        // std::cout << std::endl;
        if (local_node != UInt(-1)) {
          continue;
        }

        local_node = nodes->size();

        NodeInfo info(local_node, pos, direction);
        nodes->push_back(pos);
        nodes_global_ids->push_back(global_node);
        nodes_flags->push_back(flag | NodeFlag::_pure_ghost);
        new_nodes.push_back(info.node);
        node_list.push_back(info);

        nodes_prank[info.node] = proc;
        event.getList().push_back(local_node);
      }
    };

    /* ---------------------------------------------------------------------- */
    auto && intersections_with_right =
        bbox_left.intersection(bbox_right, *communicator);
    auto && intersections_with_left =
        bbox_right.intersection(bbox_left, *communicator);

    std::vector<CommunicationRequest> send_requests;
    std::vector<std::unique_ptr<DynamicCommunicationBuffer>> send_buffers;

    // sending nodes in the common zones
    auto send_intersections = [&](auto & intersections, auto send_count) {
      for (auto && data : intersections) {
        auto proc = std::get<0>(data);

        // Send local nodes if intersects with remote
        const auto & intersection_with_proc = std::get<1>(data);
        if (intersection_with_proc) {
          send_requests.push_back(
              extract_and_send_nodes(intersection_with_proc, nodes_right,
                                     send_buffers, proc, send_count));
        }

        send_count += 2;
      }
    };

    auto recv_intersections = [&](auto & intersections, auto recv_count) {
      for (auto && data : intersections) {
        auto proc = std::get<0>(data);

        // receive remote nodes if intersects with local
        const auto & intersection_with_proc = std::get<1>(data);
        if (intersection_with_proc) {
          recv_and_extract_nodes(nodes_right, proc, recv_count);
        }

        recv_count += 2;
      }
    };

    send_intersections(intersections_with_left, 0);
    send_intersections(intersections_with_right, 1);

    recv_intersections(intersections_with_right, 0);
    recv_intersections(intersections_with_right, 1);

    Communicator::waitAll(send_requests);
    Communicator::freeCommunicationRequest(send_requests);

    this->sendEvent(event);
  } // end distributed work

  auto to_sort = [&](auto && info1, auto && info2) -> bool {
    return info1.position < info2.position;
  };

  // sort nodes based on their distance to lower corner
  std::sort(nodes_left.begin(), nodes_left.end(), to_sort);
  std::sort(nodes_right.begin(), nodes_right.end(), to_sort);

  // function to change the master of nodes
  auto updating_master = [&](auto & old_master, auto & new_master) {
    if (old_master == new_master) {
      return;
    }

    auto slaves = periodic_master_slave.equal_range(old_master);
    AKANTU_DEBUG_ASSERT(
        isPeriodicMaster(
            old_master), // slaves.first != periodic_master_slave.end(),
        "Cannot update master " << old_master << ", its not a master node!");
    decltype(periodic_master_slave) tmp_master_slave;
    for (auto it = slaves.first; it != slaves.second; ++it) {
      auto slave = it->second;
      tmp_master_slave.insert(std::make_pair(new_master, slave));
      periodic_slave_master[slave] = new_master;
    }

    periodic_master_slave.erase(old_master);
    (*nodes_flags)[old_master] &= ~NodeFlag::_periodic_master;
    addPeriodicSlave(old_master, new_master);

    for (auto && data : tmp_master_slave) {
      addPeriodicSlave(data.second, data.first);
    }
  };

  // handling 2 nodes that are periodic
  auto match_found = [&](auto & info1, auto & info2) {
    const auto & node1 = info1.node;
    const auto & node2 = info2.node;

    auto master = node1;
    bool node1_side_master = false;
    if (isPeriodicMaster(node1)) {
      node1_side_master = true;
    } else if (isPeriodicSlave(node1)) {
      node1_side_master = true;
      master = periodic_slave_master[node1];
    }

    auto node2_master = node2;
    if (isPeriodicSlave(node2)) {
      node2_master = periodic_slave_master[node2];
    }

    if (node1_side_master) {
      if (isPeriodicSlave(node2)) {
        updating_master(node2_master, master);
        return;
      }

      if (isPeriodicMaster(node2)) {
        updating_master(node2, master);
        return;
      }

      addPeriodicSlave(node2, master);
    } else {
      if (isPeriodicSlave(node2)) {
        addPeriodicSlave(node1, node2_master);
        return;
      }

      if (isPeriodicMaster(node2)) {
        addPeriodicSlave(node1, node2);
        return;
      }

      addPeriodicSlave(node2, node1);
    }
  };

  // matching the nodes from 2 lists
  auto match_pairs = [&](auto & nodes_1, auto & nodes_2) {
    // Guillaume to Nico: It seems that the list of nodes is not sorted
    // as it was: therefore the loop cannot be truncated anymore.
    // Otherwise many pairs are missing.
    // I replaced (temporarily?) for the N^2 loop so as not to miss
    // any pbc pair.
    //

    // auto it = nodes_2.begin();
    // for every nodes in 1st list
    for (auto && info1 : nodes_1) {
      auto & pos1 = info1.position;
      // auto it_cur = it;

      // try to find a match in 2nd list
      for (auto && info2 : nodes_2) {
        // auto & info2 = *it_cur;
        auto & pos2 = info2.position;

        auto dist = pos1.distance(pos2) / length;
        if (dist < tolerance) {
          // handles the found matches
          match_found(info1, info2);
          // it = it_cur;
          break;
        }
      }
    }
  };

  match_pairs(nodes_left, nodes_right);
  // match_pairs(nodes_right, nodes_left);

  this->updatePeriodicSynchronizer();

  this->is_periodic = true;
}

/* -------------------------------------------------------------------------- */
void Mesh::wipePeriodicInfo() {
  this->is_periodic = false;

  this->periodic_slave_master.clear();
  this->periodic_master_slave.clear();

  for (auto && flags : *nodes_flags) {
    flags &= ~NodeFlag::_periodic_mask;
  }
}

/* -------------------------------------------------------------------------- */
void Mesh::updatePeriodicSynchronizer() {
  if (not this->periodic_node_synchronizer) {
    this->periodic_node_synchronizer =
        std::make_unique<PeriodicNodeSynchronizer>(
            *this, this->getID() + ":periodic_synchronizer",
            this->getMemoryID(), false);
  }

  this->periodic_node_synchronizer->update();
}

} // namespace akantu
