/**
 * @file mesh_abstract_intersector.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Mon Jul 13 2015
 * @date last modification: Mon Jul 13 2015
 *
 * @brief Abstract class for intersection computations
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

#ifndef __AKANTU_MESH_ABSTRACT_INTERSECTOR_HH__
#define __AKANTU_MESH_ABSTRACT_INTERSECTOR_HH__

#include "aka_common.hh"
#include "mesh_geom_abstract.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/**
 * @brief Class used to perform intersections on a mesh and construct output data
 */
template<class Query>
class MeshAbstractIntersector : public MeshGeomAbstract {

public:
  /// Construct from mesh
  explicit MeshAbstractIntersector(Mesh & mesh);

  /// Destructor
  virtual ~MeshAbstractIntersector();

public:
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

  /// Compute whatever result is needed from the user
  virtual void buildResultFromQueryList(const std::list<Query> & query_list) = 0;
};

__END_AKANTU__

#include "mesh_abstract_intersector_tmpl.hh"

#endif // __AKANTU_MESH_ABSTRACT_INTERSECTOR_HH__