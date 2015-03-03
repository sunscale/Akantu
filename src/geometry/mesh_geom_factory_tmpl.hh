/**
 * @file   mesh_geom_factory_tmpl.hh
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Thu Feb 26 2015
 * @date last modification: Mon Mar 2 2015
 *
 * @brief  Class for constructing the CGAL primitives of a mesh
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

#ifndef _AKANTU_MESH_GEOM_FACTORY_TMPL_HH__
#define _AKANTU_MESH_GEOM_FACTORY_TMPL_HH__

#include "mesh_geom_factory.hh"

#include "aka_common.hh"
#include "tree_type_helper.hh"
#include "triangle.hh"

#include <CGAL/Cartesian.h>

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

typedef CGAL::Cartesian<Real> K;

template<UInt d, ElementType type>
MeshGeomFactory<d, type>::MeshGeomFactory(const Mesh & mesh) :
  MeshGeomAbstract(mesh),
  data_tree(NULL),
  primitive_list()
{}

template<UInt d, ElementType type>
MeshGeomFactory<d, type>::~MeshGeomFactory() {
  delete data_tree;
}

template<>
void MeshGeomFactory<2, _triangle_3>::constructData() {
  /// Variable declaration instead of template arguments
  const UInt d = 2;
  const ElementType type = _triangle_3;
  const GhostType ghost_type = _not_ghost;

  primitive_list.clear();

  UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);
  const Array<UInt> & connectivity = mesh.getConnectivity(type, ghost_type);
  const Array<Real> & nodes = mesh.getNodes();

  Array<UInt>::const_vector_iterator begin = connectivity.begin(nb_nodes_per_element);
  Array<UInt>::const_vector_iterator it    = connectivity.begin(nb_nodes_per_element);
  Array<UInt>::const_vector_iterator end   = connectivity.end(nb_nodes_per_element);

  /// This loop builds the list of triangle primitives
  for (; it != end ; ++it) {
    const Vector<UInt> & el_connectivity = *it;
    Matrix<Real> node_coordinates(d, nb_nodes_per_element);

    for (UInt i = 0 ; i < nb_nodes_per_element ; i++)
      for (UInt j = 0 ; j < d ; j++)
        node_coordinates(j, i) = nodes(el_connectivity(i), j);

    /// There's a problem here with constructors :
    ///   - 2D constructors take 2 arguments
    ///   - 3D constructors take 3 arguments
    /// Need for a solution that can treat the two cases together
    TreeTypeHelper<d, type>::point_type a(node_coordinates(0, 0), node_coordinates(1, 0), 0.);
    TreeTypeHelper<d, type>::point_type b(node_coordinates(0, 1), node_coordinates(1, 1), 0.);
    TreeTypeHelper<d, type>::point_type c(node_coordinates(0, 2), node_coordinates(1, 2), 0.);

    /// Number of arguments here should not be a problem with polygon/polyhedra
    /// Also it is possible to create tree or lists from polyhedric surfaces
    Triangle<K> t(a, b, c);
    //TreeTypeHelper<d, type>::primitive_type t(a, b, c);
    UInt id = it - begin;
    t.setId(id);
    primitive_list.push_back(t);
  }

  /// Construct the tree
  if (data_tree) {
    delete data_tree;
    data_tree = NULL;
  }

  data_tree = new TreeTypeHelper<d, type>::tree(primitive_list.begin(), primitive_list.end());
}

template<UInt d, ElementType el_type>
UInt MeshGeomFactory<d, el_type>::numberOfIntersectionsWithInterface(const K::Segment_3 & interface) const {
  return data_tree->number_of_intersected_primitives(interface);
}

__END_AKANTU__

#endif // _AKANTU_MESH_GEOM_FACTORY_TMPL_HH__
