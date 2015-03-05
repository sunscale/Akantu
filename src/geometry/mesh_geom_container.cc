/**
 * @file   mesh_geom_container.cc
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Fri Feb 27 2015
 * @date last modification: Thu Mar 5 2015
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

#include "mesh_geom_container.hh"
#include "mesh_geom_factory.hh"
#include "mesh.hh"

#include <CGAL/Cartesian.h>

__BEGIN_AKANTU__

typedef CGAL::Cartesian<Real> K;

MeshGeomContainer::MeshGeomContainer(const Mesh & mesh):
  MeshGeomAbstract(mesh)
{}

MeshGeomContainer::~MeshGeomContainer() {
  GeomMap::type_iterator it, end;

  it  = constructor_map.firstType();
  end = constructor_map.lastType();

  for (; it != end ; it++) {
    MeshGeomAbstract * p = constructor_map(*it);

    if (p) {
      delete p;
    }
  }
}

void MeshGeomContainer::constructData() {
  const UInt spatial_dim = mesh.getSpatialDimension();

  Mesh::type_iterator it = mesh.firstType(spatial_dim, _not_ghost);
  Mesh::type_iterator end = mesh.lastType(spatial_dim, _not_ghost);

  /// Loop over the element types of the mesh and construct the primitive trees
  for (; it != end ; ++it) {
    switch(spatial_dim) {
      case 1:
        AKANTU_DEBUG_WARNING("Geometry in 1D is undefined");
        break;

      case 2:
        switch(*it) {
          case _triangle_3: {
            MeshGeomFactory<2, _triangle_3> * factory = new MeshGeomFactory<2, _triangle_3>(mesh);
            constructor_map(factory, _triangle_3);
            factory->constructData();
            break;
          }

          /// Need to implement these cases when needed
          case _not_defined:
          case _point_1:
          case _segment_2:
          case _segment_3:
          case _triangle_6: // 2nd order element : harder to implement (cf Gomes & Awruch 2001)
          case _tetrahedron_4: // needs implementing
          case _tetrahedron_10:
          case _quadrangle_4: // may need implementing
          case _quadrangle_8: // 2nd order element : harder to implement (cf Gomes & Awruch 2001)
          case _hexahedron_8: // may need implementing
          case _pentahedron_6:
          case _max_element_type:
            break;
        }
        break;
      
      case 3:
        break;
    }
  }
}

UInt MeshGeomContainer::numberOfIntersectionsWithInterface(const K::Segment_3 & interface) const {
  UInt total = 0;

  GeomMap::type_iterator it = constructor_map.firstType();
  GeomMap::type_iterator end = constructor_map.lastType();

  for (; it != end ; ++it) {
    total += constructor_map(*it)->numberOfIntersectionsWithInterface(interface);
  }

  return total;
}

Mesh * MeshGeomContainer::computeIntersectionWithLinearInterface(const K::Segment_3 & interface) {
  return NULL;
}

const MeshGeomAbstract * MeshGeomContainer::getFactoryForElementType(ElementType el_type) const {
  return constructor_map(el_type);
}

__END_AKANTU__
