/**
 * @file mesh_segment_intersector.hh
 *
 * @author Lucas Frerot
 *
 * @date creation: Wed Apr 29 2015
 * @date last modification: Wed Apr 29 2015
 *
 * @brief Computation of mesh intersection with segments
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2015 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef __AKANTU_MESH_SEGMENT_INTERSECTOR_HH__
#define __AKANTU_MESH_SEGMENT_INTERSECTOR_HH__

#include "aka_common.hh"
#include "mesh_geom_intersector.hh"

#include <CGAL/Cartesian.h>

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/// Here, we know what kernel we have to use
typedef CGAL::Cartesian<Real> K;

template<UInt dim, ElementType type>
class MeshSegmentIntersector : public MeshGeomIntersector<dim, type, Triangle<K>, K::Segment_3, K> {
  /// Parent class type
  typedef MeshGeomIntersector<dim, type, Triangle<K>, K::Segment_3, K> parent_type;

  /// Result of intersection type
  typedef typename parent_type::result_type result_type;

  /// Pair of segments and element id
  typedef std::pair<K::Segment_3, UInt> pair_type;

public:
  /// Construct from mesh
  explicit MeshSegmentIntersector(const Mesh & mesh, Mesh & result_mesh);

  /// Destructor
  virtual ~MeshSegmentIntersector();

public:
  /**
   * @brief Computes the intersection of the mesh with a segment
   *
   * @param query the segment to compute the intersections with the mesh
   */
  virtual void computeIntersectionQuery(const K::Segment_3 & query);

  /// Compute the list of queries
  virtual void computeIntersectionQueryList(const std::list<K::Segment_3> & query_list);
  
  /// Compute the list of queries and adds sets the physical name of the elements created
  virtual void computeIntersectionQueryList(const std::list<K::Segment_3> & query_list,
                                            const std::string & physical_name);

protected:
  /// Compute segments from intersection list
  void computeSegments(const std::list<result_type> & intersections,
                       std::list<pair_type> & segments,
                       const K::Segment_3 & query);

protected:
  /// Result mesh
  Mesh & result_mesh;

  /// Physical name of the current batch of queries
  std::string current_physical_name;
};
 
__END_AKANTU__

#include "mesh_segment_intersector_tmpl.hh"

#endif // __AKANTU_MESH_SEGMENT_INTERSECTOR_HH__
