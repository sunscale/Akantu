/**
 * @file mesh_sphere_intersector.hh
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
  explicit MeshSphereIntersector(Mesh & mesh);

  /// Destructor
  virtual ~MeshSphereIntersector();

public:
  /// Construct the primitive tree object
  virtual void constructData(GhostType ghost_type = _not_ghost);

  /**
   * @brief Computes the intersection of the mesh with a sphere
   *
   * @param query (sphere) to compute the intersections with the mesh
   */
  virtual void computeIntersectionQuery(const SK::Sphere_3 & query){
    AKANTU_DEBUG_ERROR("This function is not implemented for spheres (It was to generic and has been replaced by computeMeshQueryIntersectionPoint");
  }

  /// Compute intersection points between the mesh and a query
  virtual void computeMeshQueryIntersectionPoint(const SK::Sphere_3 & query, UInt nb_old_nodes);

  /// Build the IGFEM mesh
  virtual void buildResultFromQueryList(const std::list<SK::Sphere_3> & query){
    AKANTU_DEBUG_ERROR("This function is no longer implemented to split geometrical operations and dedicated result construction");
  }

  /// Set the tolerance
  void setToleranceIntersectionOnNode(UInt tol) {
    this->tol_intersection_on_node = tol;
  }

protected:
  /// tolerance for which the intersection is considered on the mesh node (relative to the segment lenght)
  Real tol_intersection_on_node;

  /// number of primitive in an element of the template type
  UInt nb_prim_by_el;

};
 
__END_AKANTU__

#include "mesh_sphere_intersector_tmpl.hh"

#endif // __AKANTU_MESH_SPHERE_INTERSECTOR_HH__
