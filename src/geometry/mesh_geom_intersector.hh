/**
 * @file mesh_geom_intersector.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Wed Apr 29 2015
 * @date last modification: Wed Apr 29 2015
 *
 * @brief General class for intersection computations
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

#ifndef __AKANTU_MESH_GEOM_INTERSECTOR_HH__
#define __AKANTU_MESH_GEOM_INTERSECTOR_HH__

#include "aka_common.hh"
#include "mesh_geom_abstract.hh"
#include "mesh_geom_factory.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/**
 * @brief Class used to perform intersections on a mesh and construct output data
 */
template<UInt dim, ElementType type, class Primitive, class Query, class Kernel>
class MeshGeomIntersector : public MeshGeomAbstract {

protected:
  /// Result of intersection function type
  typedef typename IntersectionTypeHelper<TreeTypeHelper<Primitive, Kernel>, Query>::intersection_type result_type;

public:
  /// Construct from mesh
  explicit MeshGeomIntersector(const Mesh & mesh);

  /// Destructor
  virtual ~MeshGeomIntersector();

public:
  /// Construct the primitive tree object
  virtual void constructData();

  /**
   * @brief Compute the intersection with a query object
   *
   * This function needs to be implemented for every subclass. It computes the intersections
   * with the tree of primitives and creates the data for the user.
   * 
   * @param query the CGAL primitive of the query object
   */
  virtual void computeIntersectionQuery(const Query & query) = 0;

  /// Compute list of queries
  virtual void computeIntersectionQueryList(const std::list<Query> & query_list);

protected:
  /// Factory object containing the primitive tree
  MeshGeomFactory<dim, type, Primitive, Kernel> factory;
};

__END_AKANTU__

#include "mesh_geom_intersector_tmpl.hh"

#endif // __AKANTU_MESH_GEOM_INTERSECTOR_HH__
