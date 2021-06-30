/**
 * @file   mesh_utils_pbc.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Wed Feb 09 2011
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  periodic boundary condition connectivity tweak
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
#include <map>
/* -------------------------------------------------------------------------- */
#include "element_group.hh"
#include "mesh_accessor.hh"
#include "mesh_utils.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
/// class that sorts a set of nodes of same coordinates in 'dir' direction
class CoordinatesComparison {
public:
  CoordinatesComparison(const UInt dimension, const UInt dir_1,
                        const UInt dir_2, Real normalization, Real tolerance,
                        const Array<Real> & coords)
      : dim(dimension), dir_1(dir_1), dir_2(dir_2),
        normalization(normalization), tolerance(tolerance),
        coords_it(coords.begin(dim)) {}
  // answers the question whether n2 is larger or equal to n1
  bool operator()(const UInt n1, const UInt n2) {
    Vector<Real> coords_n1 = coords_it[n1];
    Vector<Real> coords_n2 = coords_it[n2];
    return this->operator()(coords_n1, coords_n2);
  }

  bool operator()(const Vector<Real> & coords_n1,
                  const Vector<Real> & coords_n2) const {
    Real diff = coords_n1(dir_1) - coords_n2(dir_1);
    ;
    if (dim == 2 || std::abs(diff) / normalization > tolerance) {
      return (diff <= 0.);
    }
    if (dim > 2) {
      diff = coords_n1(dir_2) - coords_n2(dir_2);
      return (diff <= 0);
    }
    return true;
  }

private:
  UInt dim;
  UInt dir_1;
  UInt dir_2;
  Real normalization;
  Real tolerance;
  const Array<Real>::const_vector_iterator coords_it;
};

/* -------------------------------------------------------------------------- */
void MeshUtils::computePBCMap(const Mesh & mesh, const UInt dir,
                              std::map<UInt, UInt> & pbc_pair) {
  Array<UInt> selected_left;
  Array<UInt> selected_right;

  const UInt dim = mesh.getSpatialDimension();
  auto it = mesh.getNodes().begin(dim);
  auto end = mesh.getNodes().end(dim);

  if (dim <= dir) {
    return;
  }

  const Vector<Real> & lower_bounds = mesh.getLowerBounds();
  const Vector<Real> & upper_bounds = mesh.getUpperBounds();

  AKANTU_DEBUG_INFO("min " << lower_bounds(dir));
  AKANTU_DEBUG_INFO("max " << upper_bounds(dir));

  for (UInt node = 0; it != end; ++it, ++node) {
    const Vector<Real> & coords = *it;
    AKANTU_DEBUG_TRACE("treating " << coords(dir));
    if (Math::are_float_equal(coords(dir), lower_bounds(dir))) {
      AKANTU_DEBUG_TRACE("pushing node " << node << " on the left side");
      selected_left.push_back(node);
    } else if (Math::are_float_equal(coords(dir), upper_bounds(dir))) {
      selected_right.push_back(node);
      AKANTU_DEBUG_TRACE("pushing node " << node << " on the right side");
    }
  }

  AKANTU_DEBUG_INFO("found "
                    << selected_left.size() << " and " << selected_right.size()
                    << " nodes at each boundary for direction " << dir);

  // match pairs
  MeshUtils::matchPBCPairs(mesh, dir, selected_left, selected_right, pbc_pair);
}

/* -------------------------------------------------------------------------- */
void MeshUtils::computePBCMap(const Mesh & mesh,
                              const std::pair<ID, ID> & surface_pair,
                              std::map<UInt, UInt> & pbc_pair) {

  Array<UInt> selected_first;
  Array<UInt> selected_second;

  // find nodes on surfaces
  const ElementGroup & first_surf = mesh.getElementGroup(surface_pair.first);
  const ElementGroup & second_surf = mesh.getElementGroup(surface_pair.second);

  // if this surface pair is not on this proc
  if (first_surf.getNbNodes() == 0 || second_surf.getNbNodes() == 0) {
    AKANTU_DEBUG_WARNING("computePBCMap has at least one surface without any "
                         "nodes. I will ignore it.");
    return;
  }

  // copy nodes from element group
  selected_first.copy(first_surf.getNodeGroup().getNodes());
  selected_second.copy(second_surf.getNodeGroup().getNodes());

  // coordinates
  const Array<Real> & coords = mesh.getNodes();
  const UInt dim = mesh.getSpatialDimension();

  // variables to find min and max of surfaces
  Real first_max[3];
  Real first_min[3];
  Real second_max[3];
  Real second_min[3];
  for (UInt i = 0; i < dim; ++i) {
    first_min[i] = std::numeric_limits<Real>::max();
    second_min[i] = std::numeric_limits<Real>::max();
    first_max[i] = -std::numeric_limits<Real>::max();
    second_max[i] = -std::numeric_limits<Real>::max();
  }

  // find min and max of surface nodes
  for (auto it = selected_first.begin(); it != selected_first.end(); ++it) {
    for (UInt i = 0; i < dim; ++i) {
      if (first_min[i] > coords(*it, i)) {
        first_min[i] = coords(*it, i);
      }
      if (first_max[i] < coords(*it, i)) {
        first_max[i] = coords(*it, i);
      }
    }
  }
  for (auto it = selected_second.begin(); it != selected_second.end(); ++it) {
    for (UInt i = 0; i < dim; ++i) {
      if (second_min[i] > coords(*it, i)) {
        second_min[i] = coords(*it, i);
      }
      if (second_max[i] < coords(*it, i)) {
        second_max[i] = coords(*it, i);
      }
    }
  }

  // find direction of pbc
  Int first_dir = -1;
#ifndef AKANTU_NDEBUG
  Int second_dir = -2;
#endif
  for (UInt i = 0; i < dim; ++i) {
    if (Math::are_float_equal(first_min[i], first_max[i])) {
      first_dir = i;
    }
#ifndef AKANTU_NDEBUG
    if (Math::are_float_equal(second_min[i], second_max[i])) {
      second_dir = i;
    }
#endif
  }

  AKANTU_DEBUG_ASSERT(first_dir == second_dir,
                      "Surface pair has not same direction. Surface "
                          << surface_pair.first << " dir=" << first_dir
                          << " ; Surface " << surface_pair.second
                          << " dir=" << second_dir);
  UInt dir = first_dir;

  // match pairs
  if (first_min[dir] < second_min[dir]) {
    MeshUtils::matchPBCPairs(mesh, dir, selected_first, selected_second,
                             pbc_pair);
  } else {
    MeshUtils::matchPBCPairs(mesh, dir, selected_second, selected_first,
                             pbc_pair);
  }
}

/* -------------------------------------------------------------------------- */
void MeshUtils::matchPBCPairs(const Mesh & mesh, const UInt dir,
                              Array<UInt> & selected_left,
                              Array<UInt> & selected_right,
                              std::map<UInt, UInt> & pbc_pair) {

  // tolerance is that large because most meshers generate points coordinates
  // with single precision only (it is the case of GMSH for instance)
  Real tolerance = 1e-7;
  const UInt dim = mesh.getSpatialDimension();
  Real normalization = mesh.getUpperBounds()(dir) - mesh.getLowerBounds()(dir);

  AKANTU_DEBUG_ASSERT(std::abs(normalization) > Math::getTolerance(),
                      "In matchPBCPairs: The normalization is zero. "
                          << "Did you compute the bounding box of the mesh?");

  auto odir_1 = UInt(-1);
  auto odir_2 = UInt(-1);

  if (dim == 3) {
    if (dir == _x) {
      odir_1 = _y;
      odir_2 = _z;
    } else if (dir == _y) {
      odir_1 = _x;
      odir_2 = _z;
    } else if (dir == _z) {
      odir_1 = _x;
      odir_2 = _y;
    }
  } else if (dim == 2) {
    if (dir == _x) {
      odir_1 = _y;
    } else if (dir == _y) {
      odir_1 = _x;
    }
  }

  CoordinatesComparison compare_nodes(dim, odir_1, odir_2, normalization,
                                      tolerance, mesh.getNodes());

  std::sort(selected_left.begin(), selected_left.end(), compare_nodes);
  std::sort(selected_right.begin(), selected_right.end(), compare_nodes);

  auto it_left = selected_left.begin();
  auto end_left = selected_left.end();

  auto it_right = selected_right.begin();
  auto end_right = selected_right.end();

  auto nit = mesh.getNodes().begin(dim);

  while ((it_left != end_left) && (it_right != end_right)) {
    UInt i1 = *it_left;
    UInt i2 = *it_right;

    Vector<Real> coords1 = nit[i1];
    Vector<Real> coords2 = nit[i2];

    AKANTU_DEBUG_TRACE("do I pair? " << i1 << "(" << coords1 << ") with" << i2
                                     << "(" << coords2 << ") in direction "
                                     << dir);

    Real dx = 0.0;
    Real dy = 0.0;
    if (dim >= 2) {
      dx = coords1(odir_1) - coords2(odir_1);
    }
    if (dim == 3) {
      dy = coords1(odir_2) - coords2(odir_2);
    }

    if (std::abs(dx * dx + dy * dy) / normalization < tolerance) {
      // then i match these pairs
      if (pbc_pair.count(i2) != 0U) {
        i2 = pbc_pair[i2];
      }
      pbc_pair[i1] = i2;

      AKANTU_DEBUG_TRACE("pairing " << i1 << "(" << coords1 << ") with" << i2
                                    << "(" << coords2 << ") in direction "
                                    << dir);
      ++it_left;
      ++it_right;
    } else if (compare_nodes(coords1, coords2)) {
      ++it_left;
    } else {
      ++it_right;
    }
  }
  AKANTU_DEBUG_INFO("found " << pbc_pair.size() << " pairs for direction "
                             << dir);
}

} // namespace akantu
