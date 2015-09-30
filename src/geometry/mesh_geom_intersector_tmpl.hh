/**
 * @file mesh_geom_intersector_tmpl.hh
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

#ifndef __AKANTU_MESH_GEOM_INTERSECTOR_TMPL_HH__
#define __AKANTU_MESH_GEOM_INTERSECTOR_TMPL_HH__

#include "aka_common.hh"
#include "mesh_geom_intersector.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

template<UInt dim, ElementType type, class Primitive, class Query, class Kernel>
MeshGeomIntersector<dim, type, Primitive, Query, Kernel>::MeshGeomIntersector(Mesh & mesh,
									      const ID & id,
									      const MemoryID & memory_id) :
									      MeshAbstractIntersector<Query>(mesh, id, memory_id),
  factory(mesh)
{}

template<UInt dim, ElementType type, class Primitive, class Query, class Kernel>
MeshGeomIntersector<dim, type, Primitive, Query, Kernel>::~MeshGeomIntersector()
{}

template<UInt dim, ElementType type, class Primitive, class Query, class Kernel>
void MeshGeomIntersector<dim, type, Primitive, Query, Kernel>::constructData() {
  this->new_nodes->resize(0);
  factory.constructData();
}

__END_AKANTU__

#endif // __AKANTU_MESH_GEOM_INTERSECTOR_TMPL_HH__

