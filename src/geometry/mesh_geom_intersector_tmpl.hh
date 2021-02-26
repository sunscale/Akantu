/**
 * @file   mesh_geom_intersector_tmpl.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
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

#ifndef AKANTU_MESH_GEOM_INTERSECTOR_TMPL_HH_
#define AKANTU_MESH_GEOM_INTERSECTOR_TMPL_HH_

#include "aka_common.hh"
#include "mesh_geom_intersector.hh"

/* -------------------------------------------------------------------------- */

namespace akantu {

template <UInt dim, ElementType type, class Primitive, class Query,
          class Kernel>
MeshGeomIntersector<dim, type, Primitive, Query, Kernel>::MeshGeomIntersector(
    Mesh & mesh)
    : MeshAbstractIntersector<Query>(mesh), factory(mesh) {}

template <UInt dim, ElementType type, class Primitive, class Query,
          class Kernel>
void MeshGeomIntersector<dim, type, Primitive, Query, Kernel>::constructData(
    GhostType ghost_type) {
  this->intersection_points->resize(0);
  factory.constructData(ghost_type);
}

} // namespace akantu

#endif // AKANTU_MESH_GEOM_INTERSECTOR_TMPL_HH_
