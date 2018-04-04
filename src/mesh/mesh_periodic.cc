/**
 * @file   mesh_pbc.cc
 *
 * @author Nicolas Richart
 *
 * @date creation  Sat Feb 10 2018
 *
 * @brief Implementation of the perdiodicity capabilities in the mesh
 *
 * @section LICENSE
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
#include "mesh.hh"
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

  std::cout << bbox << std::endl;

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

namespace {
  struct NodeInfo {
    NodeInfo() {}
    NodeInfo(UInt spatial_dimension) : position(spatial_dimension) {}
    NodeInfo(UInt node, const Vector<Real> & position,
             const SpatialDirection & direction)
        : node(node), position(position) {
      this->direction_position = position(direction);
      this->position(direction) = 0.;
    }

    NodeInfo(const NodeInfo & other)
        : node(other.node), position(other.position),
          direction_position(other.direction_position) {}

    UInt node{0};
    Vector<Real> position;
    Real direction_position{0.};
  };

  // std::ostream & operator<<(std::ostream & stream, const NodeInfo & info) {
  //   stream << info.node << " " << info.position << " "
  //          << info.direction_position;
  //   return stream;
  // }
}

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
  std::cout << prank << " - left:" << list_left.size()
            << " - right:" << list_right.size() << std::endl;

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
    return info;
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
    auto extract_and_send_nodes = [&](const auto & bbox, const auto & node_list,
                                      auto & send_buffers, auto proc,
                                      auto cnt) {
      send_buffers.emplace_back();
      auto & buffer = send_buffers.back();

      for (auto & info : node_list) {
        if (bbox.contains(info.position) and isLocalOrMasterNode(info.node)) {
          Vector<Real> pos = info.position;
          pos(direction) = info.direction_position;
          buffer << getNodeGlobalId(info.node);
          buffer << pos;
        }
      }

      auto tag = Tag::genTag(prank, cnt, Tag::_PERIODIC_NODES);
      return communicator->asyncSend(buffer, proc, tag);
    };

    auto recv_and_extract_nodes = [&](auto & bbox, auto & node_list,
                                      auto & buffer, const auto proc,
                                      auto cnt) {
      if (not bbox)
        return;

      buffer.reset();
      auto tag = Tag::genTag(proc, cnt, Tag::_PERIODIC_NODES);
      communicator->receive(buffer, proc, tag);

      while (not buffer.empty()) {
        NodeInfo info(spatial_dimension);
        Vector<Real> pos(spatial_dimension);
        UInt global_node;
        buffer >> global_node;

        buffer >> pos;
        info.position = pos;

        info.direction_position = pos(direction);
        info.position(direction) = 0;
        // info.distance = lower_bound.distance(info.position);

        info.node = getNodeLocalId(global_node);
        if (info.node != UInt(-1))
          continue;

        info.node = nodes->size();

        nodes->push_back(pos);
        nodes_global_ids->push_back(global_node);
        nodes_type.push_back(NodeType(-100));
        new_nodes.push_back(info.node);
        node_list.push_back(info);

        nodes_prank[info.node] = proc;
      }
    };

    auto && intersection_right =
        bbox_left.intersection(bbox_right, *communicator);
    auto && intersection_left =
        bbox_right.intersection(bbox_left, *communicator);

    auto prank = communicator->whoAmI();
    auto nb_proc = communicator->getNbProc();

    using buffers_t = std::vector<DynamicCommunicationBuffer>;
    std::vector<CommunicationRequest> send_requests;
    buffers_t send_buffers;

    for (auto && data :
         zip(arange(nb_proc), intersection_left, intersection_right)) {
      auto proc = std::get<0>(data);
      if (proc == prank)
        continue;

      buffers_t::iterator it;

      const auto & bbox_p_send_left = std::get<1>(data);

      if (bbox_p_send_left) {
        send_requests.push_back(extract_and_send_nodes(
            bbox_p_send_left, nodes_left, send_buffers, proc, 0));
      }

      const auto & bbox_p_send_right = std::get<2>(data);
      if (bbox_p_send_right) {
        send_requests.push_back(extract_and_send_nodes(
            bbox_p_send_right, nodes_right, send_buffers, proc, 1));
      }
    }

    DynamicCommunicationBuffer buffer;
    for (auto && data :
         zip(arange(nb_proc), intersection_left, intersection_right)) {
      auto proc = std::get<0>(data);
      if (proc == prank)
        continue;

      const auto & bbox_p_send_left = std::get<1>(data);
      recv_and_extract_nodes(bbox_p_send_left, nodes_left, buffer, proc, 0);

      const auto & bbox_p_send_right = std::get<2>(data);
      recv_and_extract_nodes(bbox_p_send_right, nodes_right, buffer, proc, 1);
    }

    communicator->waitAll(send_requests);
    communicator->freeCommunicationRequest(send_requests);
  }

  auto to_sort = [&](auto && info1, auto && info2) -> bool {
    return info1.position < info2.position;
  };

  std::sort(nodes_left.begin(), nodes_left.end(), to_sort);
  std::sort(nodes_right.begin(), nodes_right.end(), to_sort);

  auto match_found = [&](auto & info1, auto & info2) {
    auto node1 = info1.node;
    auto node2 = info2.node;
    if (info1.direction_position < info2.direction_position) {
      std::swap(node1, node2);
    }

    auto mits = periodic_master_slave.equal_range(node2);
    if (mits.first != periodic_master_slave.end()) {
      for (auto mit = mits.first; mit != mits.second; ++mit) {
        periodic_master_slave.insert(std::make_pair(node1, mit->second));
        // \TODO tell processor prank[mit->second] that master is now node1
        // instead of node2
        periodic_slave_master[mit->second] = node1;
        // \TODO tell processor prank[node1] that it also has a slave on
        // prank[mit->second]
      }
      periodic_master_slave.erase(node2);
    }

    auto node1_slaves = periodic_master_slave.equal_range(node1);
    auto slave_it =
        std::find_if(node1_slaves.first, node1_slaves.second,
                     [&](auto & pair) { return pair.second == node2; });
    if (slave_it == node1_slaves.second)
      periodic_master_slave.insert(std::make_pair(node1, node2));

    periodic_slave_master[node2] = node1;
  };

  auto match_pairs = [&](auto & nodes_1, auto & nodes_2) {
    auto it = nodes_2.begin();
    for (auto && info1 : nodes_1) {
      if (not isLocalOrMasterNode(info1.node))
        continue;

      auto & pos1 = info1.position;
      auto it_cur = it;

      for (; it_cur != nodes_2.end(); ++it_cur) {
        auto & info2 = *it_cur;
        auto & pos2 = info2.position;
        auto dist = pos1.distance(pos2) / length;

        if (dist < tolerance) {
          match_found(info1, *it_cur);
          it = it_cur;
          break;
        }
      }
    }
  };

  match_pairs(nodes_left, nodes_right);
  match_pairs(nodes_right, nodes_left);

  std::cout << periodic_slave_master.size() << std::endl;

  for (auto & pair : periodic_master_slave) {
    std::cout << prank << " - " << pair.first << " " << pair.second
              << std::endl;
  }

  // std::cout << prank << " - Left" << std::endl;
  // for (auto && data : nodes_left) {
  //   std::cout << prank << " - " << data << " -- " << getNodeType(data.node)
  //             << std::endl;
  // }

  // std::cout << prank << " - Right" << std::endl;
  // for (auto && data : nodes_right) {
  //   std::cout << prank << " - " << data << " -- " << getNodeType(data.node)
  //             << std::endl;
  //}

  this->is_periodic |= 1 < direction;
}

} // akantu
