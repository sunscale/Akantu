/**
 * @file mesh_segment_intersector.hh
 *
 * @author Clement Roux-Langlois <clement.roux@epfl.ch>
 *
 * @date creation: Wed Jun 10 2015
 *
 * @brief Computation of mesh intersection with sphere(s)
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

#ifndef __AKANTU_MESH_SPHERE_INTERSECTOR_HH__
#define __AKANTU_MESH_SPHERE_INTERSECTOR_HH__

#include "aka_common.hh"
#include "mesh_geom_intersector.hh"

#include "mesh_geom_common.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/// Here, we know what kernel we have to use
typedef Spherical SK;

template<UInt dim, ElementType type>
class MeshSphereIntersector : public MeshGeomIntersector<dim, type, Line_arc<SK>, SK::Sphere_3, SK> {
  /// Parent class type
  typedef MeshGeomIntersector<dim, type, Line_arc<SK>, SK::Sphere_3, SK> parent_type;

  /// Result of intersection function type
  typedef typename IntersectionTypeHelper<TreeTypeHelper< Triangle<K>, K>, K::Segment_3>::intersection_type result_type;

  /// Pair of intersection points and element id
  typedef std::pair<SK::Circular_arc_point_3, UInt> pair_type;

public:
  /// Construct from mesh
  explicit MeshSphereIntersector(const Mesh & mesh);

  /// Destructor
  virtual ~MeshSphereIntersector();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the new_node_per_elem array
  AKANTU_GET_MACRO(NewNodePerElem, new_node_per_elem, const Array<UInt>)

public:
  /**
   * @brief Computes the intersection of the mesh with a segment
   *
   * @param query the segment to compute the intersections with the mesh
   */
  virtual void computeIntersectionQuery(const SK::Sphere_3 & query);

  /// Compute the list of queries
  virtual void computeIntersectionQueryList(const std::list<SK::Sphere_3> & query_list);
  
#if defined(AKANTU_IGFEM)
  /// Build the IGFEM mesh
  void buildIgfemMesh(const std::list<SK::Sphere_3> & query_list);
#endif

protected:
  /// new node per element TODO convert to ElementTypeMapArray<UInt>
  Array<UInt> new_node_per_elem;

  /// new node per element TODO convert to ElementTypeMapArray<UInt>
  const UInt nb_nodes_init;
};
 
__END_AKANTU__

#include "mesh_sphere_intersector_tmpl.hh"

#endif // __AKANTU_MESH_SPHERE_INTERSECTOR_HH__
