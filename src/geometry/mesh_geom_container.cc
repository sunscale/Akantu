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

#include "mesh_geom_container.hh"
#include "mesh_geom_factory.hh"
#include "mesh.hh"

#include <CGAL/Cartesian.h>

__BEGIN_AKANTU__

typedef CGAL::Cartesian<Real> K;

MeshGeomContainer::MeshGeomContainer(const Mesh & mesh):
  MeshGeomAbstract(mesh),
  factory_map(),
  interface_mesh()
{}

MeshGeomContainer::~MeshGeomContainer() {
  /*GeomMap::type_iterator it, end;

  it  = factory_map.firstType();
  end = factory_map.lastType();

  for (; it != end ; it++) {
    MeshGeomAbstract * p = factory_map(*it);
    delete p;
    p = NULL;
  }*/

  delete interface_mesh;
}

void MeshGeomContainer::constructData() {
  AKANTU_DEBUG_IN();

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
            factory_map(factory, _triangle_3);
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

Mesh * MeshGeomContainer::meshOfLinearInterface(const std::pair<K::Segment_3, std::string> & pair) {
  AKANTU_DEBUG_IN();

  std::list<Mesh *> meshes;

  GeomMap::type_iterator it = factory_map.firstType();
  GeomMap::type_iterator end = factory_map.lastType();

  for (; it != end ; ++it) {
    meshes.push_back(factory_map(*it)->meshOfLinearInterface(pair));
  }

  AKANTU_DEBUG_OUT();

  if (meshes.size() == 1) {
    return meshes.front();
  } else {
    return mergeMeshes(meshes);
  }
}

Mesh * MeshGeomContainer::mergeMeshes(const std::list<Mesh *> & meshes) {
  AKANTU_DEBUG_IN();

  UInt d = mesh.getSpatialDimension();

  // XXX there's gonna be a problem here with several element types in meshes (ID)
  Mesh * mergedMesh = new Mesh(d, "mergedMesh");

  Array<Real>                  mergedNodes(0, d, "mergedNodes");
  ElementTypeMapArray<UInt>    mergedConnectivities("mergedConnectivities");
  ElementTypeMapArray<UInt>    mergedAssociatedIds("mergedAssociatedIds");
  ElementTypeMapArray<Element> mergedAssociatedTypes("mergedAssociatedTypes");

  std::list<Mesh *>::const_iterator mesh_it = meshes.begin();
  std::list<Mesh *>::const_iterator mesh_end = meshes.end();

  for (; mesh_it != mesh_end ; ++mesh_it) {
    const Mesh * currentMesh = *mesh_it;

    Mesh::type_iterator type_it = currentMesh->firstType();
    Mesh::type_iterator type_end = currentMesh->lastType();

    for (; type_it != type_end ; ++type_it) {
      mergedMesh->addConnectivityType(*type_it);

      UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(*type_it);

      // XXX very repetitive code ahead, especially for associated data

      if (!mergedConnectivities.exists(*type_it)) {
        mergedConnectivities.alloc(0, nb_nodes_per_element, *type_it, _not_ghost);
      }

      const Array<UInt> & currentConnectivity = currentMesh->getConnectivity(*type_it);
      Array<UInt>::const_vector_iterator currentConn_it = currentConnectivity.begin(nb_nodes_per_element);
      Array<UInt>::const_vector_iterator currentConn_end = currentConnectivity.end(nb_nodes_per_element);

      /// Offset for node id
      Vector<UInt> offset(nb_nodes_per_element, mergedNodes.getSize());

      /// Append connectivity
      for (; currentConn_it != currentConn_end ; ++currentConn_it) {
        mergedConnectivities(*type_it).push_back(*currentConn_it + offset);
      }

      // ------------------

      if (!mergedAssociatedIds.exists(*type_it)) {
        mergedAssociatedIds.alloc(0, 1, *type_it, _not_ghost);
      }

      const Array<UInt> & currentAssociatedIds = currentMesh->getData<UInt>("associated_id", *type_it);
      Array<UInt>::const_scalar_iterator currentAssoId_it = currentAssociatedIds.begin();
      Array<UInt>::const_scalar_iterator currentAssoId_end = currentAssociatedIds.end();

      /// Append associated ids
      for (; currentAssoId_it != currentAssoId_end ; ++currentAssoId_it) {
        mergedAssociatedIds(*type_it).push_back(*currentAssoId_it);
      }

      // ------------------

      if (!mergedAssociatedTypes.exists(*type_it)) {
        mergedAssociatedTypes.alloc(0, 1, *type_it, _not_ghost);
      }

      const Array<Element> & currentAssociatedTypes = currentMesh->getData<Element>("associated_type", *type_it);
      Array<Element>::const_scalar_iterator currentAssoType_it = currentAssociatedTypes.begin();
      Array<Element>::const_scalar_iterator currentAssoType_end = currentAssociatedTypes.end();

      /// Append associated types
      for (; currentAssoType_it != currentAssoType_end ; ++currentAssoType_it) {
        mergedAssociatedTypes(*type_it).push_back(*currentAssoType_it);
      }
    }

    Array<Real>::const_vector_iterator currentNodes_it = currentMesh->getNodes().begin(d);
    Array<Real>::const_vector_iterator currentNodes_end = currentMesh->getNodes().end(d);

    /// Append nodes
    for (; currentNodes_it != currentNodes_end ; ++currentNodes_it) {
      mergedNodes.push_back(*currentNodes_it);
    }
  }

  mergedMesh->getNodes().copy(mergedNodes);

  ElementTypeMapArray<UInt>::type_iterator type_it = mergedConnectivities.firstType();
  ElementTypeMapArray<UInt>::type_iterator type_end = mergedConnectivities.lastType();

  for (; type_it != type_end ; ++type_it) {
    mergedMesh->getConnectivity(*type_it).copy(mergedConnectivities(*type_it));

    mergedMesh->registerData<UInt>("associated_id").alloc(0, 1, *type_it, _not_ghost);
    mergedMesh->registerData<Element>("associated_type").alloc(0, 1, *type_it, _not_ghost);

    mergedMesh->getData<UInt>("associated_id", *type_it).copy(mergedAssociatedIds(*type_it));
    mergedMesh->getData<Element>("associated_type", *type_it).copy(mergedAssociatedTypes(*type_it));
  }

  AKANTU_DEBUG_OUT();

  return mergedMesh;
}

Mesh & MeshGeomContainer::meshOfLinearInterfaces(const std::list<std::pair<K::Segment_3, std::string> > & interfaces) {
  if (interfaces.size() == 1) {
    interface_mesh = meshOfLinearInterface(interfaces.front());
    return *interface_mesh;
  }

  std::list<std::pair<K::Segment_3, std::string> >::const_iterator interfaces_it = interfaces.begin();
  std::list<std::pair<K::Segment_3, std::string> >::const_iterator interfaces_end = interfaces.end();

  std::list<Mesh *> meshes;

  for (; interfaces_it != interfaces_end ; ++interfaces_it) {
    meshes.push_back(meshOfLinearInterface(*interfaces_it));
  }

  interface_mesh = mergeMeshes(meshes);

  std::list<Mesh *>::iterator meshes_it = meshes.begin();
  std::list<Mesh *>::iterator meshes_end = meshes.end();

  for (; meshes_it != meshes_end ; ++meshes_it) {
    delete *meshes_it;
  }

  return *interface_mesh;
}

const MeshGeomAbstract * MeshGeomContainer::getFactoryForElementType(ElementType el_type) const {
  return factory_map(el_type);
}

__END_AKANTU__
