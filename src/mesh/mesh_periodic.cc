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
#include "communicator.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
void Mesh::makePeriodic(const SpatialDirection & direction) {
  const auto & lower_bound = this->getLowerBounds();
  const auto & upper_bound = this->getUpperBounds();

  Array<UInt> list_1;
  Array<UInt> list_2;

  const auto & positions = *nodes;

  for (auto && data : enumerate(make_view(positions, spatial_dimension))) {
    UInt node = std::get<0>(data);
    const auto & pos = std::get<1>(data);

    if (Math::are_float_equal(pos(direction), lower_bound(direction))) {
      list_1.push_back(node);
    }

    if (Math::are_float_equal(pos(direction), upper_bound(direction))) {
      list_2.push_back(node);
    }
  }

  this->makePeriodic(direction, list_1, list_2);
}

namespace {
  struct NodeInfo {
    NodeInfo() {}
    NodeInfo(UInt node, const Vector<Real> & position,
             const Vector<Real> & lower_node,
             const SpatialDirection & direction)
        : node(node), position(position) {
      direction_position = position(direction);
      this->position(direction) = 0.;
      this->distance = this->position.distance(lower_node);
    }

    NodeInfo(const NodeInfo & other)
        : node(other.node), distance(other.distance), position(other.position),
          direction_position(other.direction_position) {}

    UInt node{0};
    Real distance{-1.};
    Vector<Real> position;
    Real direction_position{0.};
  };

  // std::ostream & operator<<(std::ostream & stream, const NodeInfo & info) {
  //   stream << info.node << " " << info.position << " " << info.distance;
  //   return stream;
  // }
}

class BBox {
public:
  BBox(UInt spatial_dimension)
      : dim(spatial_dimension),
        lower_bounds(spatial_dimension, std::numeric_limits<Real>::max()),
        upper_bounds(spatial_dimension, -std::numeric_limits<Real>::max()) {}

  BBox(const BBox & other)
      : dim(other.dim), lower_bounds(other.lower_bounds),
        upper_bounds(other.upper_bounds) {}

  BBox & operator=(const BBox & other) {
    if (this != &other) {
      this->dim = dim;
      this->lower_bounds = other.lower_bounds;
      this->upper_bounds = other.upper_bounds;
    }
    return *this;
  }

  BBox & operator+=(const Vector<Real> & position) {
    for (auto s : arange(dim)) {
      lower_bounds(s) = std::min(lower_bounds(s), position(s));
      upper_bounds(s) = std::min(upper_bounds(s), position(s));
    }
    return *this;
  }

  const Vector<Real> & getLowerBounds() const { return lower_bounds; }
  const Vector<Real> & getUpperBounds() const { return upper_bounds; }

  Vector<Real> & getLowerBounds() { return lower_bounds; }
  Vector<Real> & getUpperBounds() { return upper_bounds; }

  void reset() {
    lower_bounds.set(std::numeric_limits<Real>::max());
    upper_bounds.set(-std::numeric_limits<Real>::max());
  }

protected:
  UInt dim;
  Vector<Real> lower_bounds;
  Vector<Real> upper_bounds;
};

/* -------------------------------------------------------------------------- */
void Mesh::makePeriodic(const SpatialDirection & direction,
                        const Array<UInt> & list_1,
                        const Array<UInt> & list_2) {
  Real tolerance = Math::getTolerance();

  const auto & positions = *nodes;
  auto lower_bound = this->getLowerBounds();
  auto upper_bound = this->getUpperBounds();
  auto length = upper_bound(direction) - lower_bound(direction);

  lower_bound(direction) = 0;
  upper_bound(direction) = 0;

  std::vector<NodeInfo> nodes_1(list_1.size());
  std::vector<NodeInfo> nodes_2(list_2.size());

  BBox bbox(spatial_dimension);
  auto to_position = [&](UInt node) {
    Vector<Real> pos(spatial_dimension);
    for (UInt s : arange(spatial_dimension)) {
      pos(s) = direction == s ? 0 : positions(node, s);
    }
    bbox += pos;
    return NodeInfo(node, pos, lower_bound, direction);
  };

  std::transform(list_1.begin(), list_1.end(), nodes_1.begin(), to_position);
  BBox bbox1 = bbox;

  bbox.reset();
  std::transform(list_2.begin(), list_2.end(), nodes_2.begin(), to_position);
  BBox bbox2 = bbox;

  if (is_distributed) {
    auto prank = communicator->whoAmI();
    auto nb_proc = communicator->getNbProc();
    Array<Real> bboxes(nb_proc, spatial_dimension * 4);
    auto * base = bboxes.storage() + prank * 4 * spatial_dimension;
    Vector<Real>(base + spatial_dimension * 0, spatial_dimension) =
        bbox1.getLowerBounds();
    Vector<Real>(base + spatial_dimension * 1, spatial_dimension) =
        bbox1.getUpperBounds();
    Vector<Real>(base + spatial_dimension * 2, spatial_dimension) =
        bbox2.getLowerBounds();
    Vector<Real>(base + spatial_dimension * 3, spatial_dimension) =
        bbox2.getUpperBounds();

    communicator->allGather(bboxes);

    for (auto p : arange(nb_proc)) {
      if (p == prank)
        continue;
    }
  }

  auto to_sort = [&](auto && info1, auto && info2) -> bool {
    return info1.distance < info2.distance;
  };

  std::sort(nodes_1.begin(), nodes_1.end(), to_sort);
  std::sort(nodes_2.begin(), nodes_2.end(), to_sort);

  auto it = nodes_2.begin();
  for (auto && info1 : nodes_1) {
    auto & pos1 = info1.position;
    auto it_cur = it;

    bool found = false;
    for (; it_cur != nodes_2.end(); ++it_cur) {
      auto & info2 = *it_cur;
      auto & pos2 = info2.position;
      auto dist = pos1.distance(pos2) / length;
      if (dist < tolerance) {
        found = true;
        it = it_cur;
        break;
      }
    }

    if (found) {
      auto node1 = info1.node;
      auto node2 = it_cur->node;
      if (info1.direction_position < it_cur->direction_position) {
        std::swap(node1, node2);
      }

      periodic_pairs.emplace(node1, std::make_pair(node2, direction));

      std::cout << "master: " << node1 << " - slave: " << node2 << std::endl;
    }
  }

  std::cout << periodic_pairs.size() << std::endl;

  this->is_periodic |= 1 < direction;
}

} // akantu
