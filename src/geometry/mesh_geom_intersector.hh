/**
 * @file   mesh_geom_intersector.hh
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

#ifndef AKANTU_MESH_GEOM_INTERSECTOR_HH_
#define AKANTU_MESH_GEOM_INTERSECTOR_HH_

#include "aka_common.hh"
#include "mesh_abstract_intersector.hh"
#include "mesh_geom_factory.hh"

/* -------------------------------------------------------------------------- */

namespace akantu {

/**
 * @brief Class used to perform intersections on a mesh and construct output
 * data
 */
template <UInt dim, ElementType type, class Primitive, class Query,
          class Kernel>
class MeshGeomIntersector : public MeshAbstractIntersector<Query> {

public:
  /// Construct from mesh
  explicit MeshGeomIntersector(Mesh & mesh);

  /// Destructor
  ~MeshGeomIntersector() override = default;

public:
  /// Construct the primitive tree object
  void constructData(GhostType ghost_type = _not_ghost) override;

protected:
  /// Factory object containing the primitive tree
  MeshGeomFactory<dim, type, Primitive, Kernel> factory;
};

} // namespace akantu

#include "mesh_geom_intersector_tmpl.hh"

#endif // AKANTU_MESH_GEOM_INTERSECTOR_HH_
