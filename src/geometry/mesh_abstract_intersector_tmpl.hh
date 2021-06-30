/**
 * @file   mesh_abstract_intersector_tmpl.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Clement Roux <clement.roux@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Apr 29 2015
 * @date last modification: Mon Jun 19 2017
 *
 * @brief  General class for intersection computations
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

#ifndef AKANTU_MESH_ABSTRACT_INTERSECTOR_TMPL_HH_
#define AKANTU_MESH_ABSTRACT_INTERSECTOR_TMPL_HH_

#include "aka_common.hh"
#include "mesh_abstract_intersector.hh"

namespace akantu {

template <class Query>
MeshAbstractIntersector<Query>::MeshAbstractIntersector(Mesh & mesh)
    : MeshGeomAbstract(mesh) {}

template <class Query>
void MeshAbstractIntersector<Query>::computeIntersectionQueryList(
    const std::list<Query> & query_list) {
  AKANTU_DEBUG_IN();

  auto query_it = query_list.begin();
  auto query_end = query_list.end();

  for (; query_it != query_end; ++query_it) {
    computeIntersectionQuery(*query_it);
  }

  AKANTU_DEBUG_OUT();
}

template <class Query>
void MeshAbstractIntersector<Query>::computeMeshQueryListIntersectionPoint(
    const std::list<Query> & query_list, UInt nb_old_nodes) {
  AKANTU_DEBUG_IN();

  auto query_it = query_list.begin();
  auto query_end = query_list.end();

  for (; query_it != query_end; ++query_it) {
    computeMeshQueryIntersectionPoint(*query_it, nb_old_nodes);
  }

  AKANTU_DEBUG_OUT();
}

} // namespace akantu

#endif // AKANTU_MESH_ABSTRACT_INTERSECTOR_TMPL_HH_
