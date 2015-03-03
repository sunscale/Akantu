/**
 * @file   mesh_geom_container.hh
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Fri Feb 27 2015
 * @date last modification: Mon Mar 2 2015
 *
 * @brief  Contains the CGAL representation of a mesh
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef __AKANTU_MESH_GEOM_CONTAINER__
#define __AKANTU_MESH_GEOM_CONTAINER__

#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_geom_abstract.hh"

__BEGIN_AKANTU__

class MeshGeomContainer : MeshGeomAbstract {

public:
  /// Construct from mesh
  explicit MeshGeomContainer(const Mesh & mesh);

  /// Destructor
  virtual ~MeshGeomContainer();

public:
  /// Constructs the geometric data from the mesh
  virtual void constructData();

  /// Get the factory object for an element type
  const MeshGeomAbstract * getFactoryForElementType(ElementType el_type) const;

protected:
  ElementTypeMap<MeshGeomAbstract *> constructor_map;
};

__END_AKANTU__

#endif // __AKANTU_MESH_GEOM_CONTAINER__
