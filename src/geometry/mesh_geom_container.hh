/**
 * @file   mesh_geom_container.hh
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Fri Feb 27 2015
 * @date last modification: Fri Mar 6 2015
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

#include <CGAL/Cartesian.h>

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

typedef CGAL::Cartesian<Real> K;

/**
 * @brief Class used to store mesh geometries for CGAL algorithms
 *
 * This class provides an interface for dialog with the MeshGeomFactory classes
 * that correspond to element types of the mesh.
 *
 * It contains a mesh object used in the intersection of a mesh with linear interfaces
 */
class MeshGeomContainer : MeshGeomAbstract {
  typedef ElementTypeMap<MeshGeomAbstract *> GeomMap;

public:
  /// Construct from mesh
  explicit MeshGeomContainer(const Mesh & mesh);

  /// Destructor
  virtual ~MeshGeomContainer();

public:
  /// Constructs the geometric data from the mesh
  virtual void constructData();

  /// Compute the number of intersections with geometric interface
  virtual UInt numberOfIntersectionsWithInterface(const K::Segment_3 & interface) const;

  /// Compute the intersection mesh with linear interface
  virtual void meshOfLinearInterface(const Interface & interface, Mesh & interface_mesh);

  /// Construct the interface mesh from several segments
  Mesh & meshOfLinearInterfaces(const std::list<Interface> & interfaces);

  /// Get the factory object for an element type
  const MeshGeomAbstract * getFactoryForElementType(ElementType el_type) const;

protected:
  GeomMap factory_map;

  Mesh interface_mesh;
};

__END_AKANTU__

#endif // __AKANTU_MESH_GEOM_CONTAINER__
