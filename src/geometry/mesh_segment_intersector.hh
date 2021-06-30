/**
 * @file   mesh_segment_intersector.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Clement Roux <clement.roux@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Apr 29 2015
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  Computation of mesh intersection with segments
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

#ifndef AKANTU_MESH_SEGMENT_INTERSECTOR_HH_
#define AKANTU_MESH_SEGMENT_INTERSECTOR_HH_

#include "aka_common.hh"
#include "mesh_geom_intersector.hh"

#include "mesh_geom_common.hh"

/* -------------------------------------------------------------------------- */

namespace akantu {

template <UInt dim, ElementType type>
class MeshSegmentIntersector
    : public MeshGeomIntersector<dim, type, Triangle<cgal::Cartesian>,
                                 cgal::Cartesian::Segment_3, cgal::Cartesian> {
  using K = cgal::Cartesian;
  /// Parent class type
  using parent_type =
      MeshGeomIntersector<dim, type, Triangle<K>, K::Segment_3, K>;

  /// Result of intersection function type
  using result_type =
      typename IntersectionTypeHelper<TreeTypeHelper<Triangle<K>, K>,
                                      K::Segment_3>::intersection_type;

  /// Pair of segments and element id
  using pair_type = std::pair<K::Segment_3, UInt>;

public:
  /// Construct from mesh
  explicit MeshSegmentIntersector(Mesh & mesh, Mesh & result_mesh);

  /// Destructor
  ~MeshSegmentIntersector() override = default;

public:
  /**
   * @brief Computes the intersection of the mesh with a segment
   *
   * @param query the segment to compute the intersections with the mesh
   */
  void computeIntersectionQuery(const K::Segment_3 & query) override;

  /// Compute intersection points between the mesh and a query
  void computeMeshQueryIntersectionPoint(const K::Segment_3 & query,
                                                 UInt nb_old_nodes) override;

  /// Compute the embedded mesh
  void
  buildResultFromQueryList(const std::list<K::Segment_3> & query_list) override;

  void setPhysicalName(const std::string & other) {
    current_physical_name = other;
  }

protected:
  /// Compute segments from intersection list
  void computeSegments(const std::list<result_type> & intersections,
                       std::set<pair_type, segmentPairsLess> & segments,
                       const K::Segment_3 & query);

protected:
  /// Result mesh
  Mesh & result_mesh;

  /// Physical name of the current batch of queries
  std::string current_physical_name;
};

} // namespace akantu

#include "mesh_segment_intersector_tmpl.hh"

#endif // AKANTU_MESH_SEGMENT_INTERSECTOR_HH_
