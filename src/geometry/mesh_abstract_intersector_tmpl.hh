/**
 * @file mesh_abstract_intersector_tmpl.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Mon Jul 13 2015
 * @date last modification: Mon Jul 13 2015
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

#ifndef __AKANTU_MESH_ABSTRACT_INTERSECTOR_TMPL_HH__
#define __AKANTU_MESH_ABSTRACT_INTERSECTOR_TMPL_HH__

#include "aka_common.hh"
#include "mesh_abstract_intersector.hh"

__BEGIN_AKANTU__

template<class Query>
MeshAbstractIntersector<Query>::MeshAbstractIntersector(Mesh & mesh,
							const ID & id,
							const MemoryID & memory_id):
  MeshGeomAbstract(mesh, id, memory_id),
  new_node_per_elem("new_node_per_elem", id),
  new_nodes(NULL),
  nb_seg_by_el(0)
{}

template<class Query>
MeshAbstractIntersector<Query>::~MeshAbstractIntersector()
{}

template<class Query>
void MeshAbstractIntersector<Query>::computeIntersectionQueryList(
  const std::list<Query> & query_list) {
  AKANTU_DEBUG_IN();
  
  typename std::list<Query>::const_iterator
    query_it = query_list.begin(),
    query_end = query_list.end();

  for (; query_it != query_end ; ++query_it) {
    computeIntersectionQuery(*query_it);
  }
  
  AKANTU_DEBUG_OUT();
}

template<class Query>
void MeshAbstractIntersector<Query>::computeMeshQueryListIntersectionPoint(
  const std::list<Query> & query_list) {
  AKANTU_DEBUG_IN();
  
  typename std::list<Query>::const_iterator
    query_it = query_list.begin(),
    query_end = query_list.end();

  for (; query_it != query_end ; ++query_it) {
    computeMeshQueryIntersectionPoint(*query_it);
  }
  
  AKANTU_DEBUG_OUT();
}

__END_AKANTU__

#endif // __AKANTU_MESH_ABSTRACT_INTERSECTOR_TMPL_HH__
