/**
 * @file   mesh_geom_container.cc
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

#include "aka_common.hh"

#include "mesh_geom_container.hh"
#include "mesh_geom_factory.hh"
#include "mesh.hh"

#include <CGAL/Cartesian.h>

/* -------------------------------------------------------------------------- */

#define MESH_GEOM_CASE(d, type)     \
  MeshGeomFactory<d, type> * factory = new MeshGeomFactory<d, type>(mesh);      \
  factory_map(factory, type);     \
  factory->constructData();

#define MESH_GEOM_CASE_2D(type) MESH_GEOM_CASE(2, type)
#define MESH_GEOM_CASE_3D(type) MESH_GEOM_CASE(3, type)

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

typedef CGAL::Cartesian<Real> K;


MeshGeomContainer::MeshGeomContainer(const Mesh & mesh):
  MeshGeomAbstract(mesh),
  factory_map(),
  interface_mesh(mesh.getSpatialDimension(), "mesh_geom")
{
  interface_mesh.addConnectivityType(_segment_2);
  interface_mesh.registerData<Element>("associated_element").alloc(0, 1, _segment_2);
}

MeshGeomContainer::~MeshGeomContainer()
{}

void MeshGeomContainer::constructData() {
  AKANTU_DEBUG_IN();

  const UInt spatial_dim = mesh.getSpatialDimension();

  Mesh::type_iterator it = mesh.firstType(spatial_dim, _not_ghost);
  Mesh::type_iterator end = mesh.lastType(spatial_dim, _not_ghost);

  /// Loop over the element types of the mesh and construct the primitive trees
  for (; it != end ; ++it) {
    ElementType type = *it; // for AKANTU_BOOST_ELEMENT_SWITCH macro
    switch(spatial_dim) {
      case 1:
        AKANTU_DEBUG_WARNING("Geometry in 1D is undefined");
        break;

      case 2:
        // Expand the list of elements when they are implemented
        AKANTU_BOOST_ELEMENT_SWITCH(MESH_GEOM_CASE_2D, (_triangle_3));
        break;
      
      case 3:
        //AKANTU_BOOST_ELEMENT_SWITCH(MESH_GEOM_CASE_3D, ());
        break;
    }
  }

  AKANTU_DEBUG_OUT();
}

UInt MeshGeomContainer::numberOfIntersectionsWithInterface(const K::Segment_3 & interface) const {
  AKANTU_DEBUG_IN();

  UInt total = 0;

  GeomMap::type_iterator it = factory_map.firstType();
  GeomMap::type_iterator end = factory_map.lastType();

  for (; it != end ; ++it) {
    total += factory_map(*it)->numberOfIntersectionsWithInterface(interface);
  }

  AKANTU_DEBUG_OUT();

  return total;
}

void MeshGeomContainer::meshOfLinearInterface(const K::Segment_3 & interface, Mesh & interface_mesh) {
  AKANTU_DEBUG_IN();

  GeomMap::type_iterator it = factory_map.firstType();
  GeomMap::type_iterator end = factory_map.lastType();

  for (; it != end ; ++it) {
    factory_map(*it)->meshOfLinearInterface(interface, interface_mesh);
  }

  AKANTU_DEBUG_OUT();
}

Mesh & MeshGeomContainer::meshOfLinearInterfaces(const std::list<K::Segment_3> & interfaces) {
  if (interfaces.size() == 1) {
    meshOfLinearInterface(interfaces.front(), interface_mesh);
    return interface_mesh;
  }

  std::list<K::Segment_3>::const_iterator interfaces_it = interfaces.begin();
  std::list<K::Segment_3>::const_iterator interfaces_end = interfaces.end();

  for (; interfaces_it != interfaces_end ; ++interfaces_it) {
    meshOfLinearInterface(*interfaces_it, interface_mesh);
  }

  return interface_mesh;
}

const MeshGeomAbstract * MeshGeomContainer::getFactoryForElementType(ElementType el_type) const {
  return factory_map(el_type);
}

__END_AKANTU__
