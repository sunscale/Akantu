/**
 * @file   mesh_geom_factory_tmpl.hh
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Thu Feb 26 2015
 * @date last modification: Thu Mar 5 2015
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
#include "geom_helper_functions.hh"

#include <CGAL/Cartesian.h>

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

typedef CGAL::Cartesian<Real> K;

template<UInt d, ElementType type>
MeshGeomFactory<d, type>::MeshGeomFactory(const Mesh & mesh) :
  MeshGeomAbstract(mesh),
  data_tree(NULL),
  primitive_list(),
  nodes(NULL),
  associated_id(NULL),
  associated_type(NULL)
{}

template<UInt d, ElementType type>
MeshGeomFactory<d, type>::~MeshGeomFactory() {
  delete data_tree;
  delete nodes;
  delete associated_id;
  delete associated_type;
}

// Need to implement this method for every case needed
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
    t.setId(it - begin);
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

template<UInt d, ElementType el_type>
Mesh * MeshGeomFactory<d, el_type>::computeIntersectionWithLinearInterface(const K::Segment_3 & interface) {
  UInt number_of_intersections = this->numberOfIntersectionsWithInterface(interface);

  if (!number_of_intersections)
    return NULL;

  std::list<typename TreeTypeHelper<d, el_type>::linear_intersection> list_of_intersections;
  std::list< std::pair<K::Segment_3, UInt> > list_of_segments; // Contains no duplicate elements

  data_tree->all_intersections(interface, std::back_inserter(list_of_intersections));
  this->constructSegments(list_of_intersections, list_of_segments);

  if (nodes) {
    delete nodes;
    nodes = NULL;
  }

  nodes = new Array<Real>(2 * list_of_segments.size(), d);
  Array<UInt> connectivity(list_of_segments.size(), 2);

  std::list< std::pair<K::Segment_3, UInt> >::iterator it  = list_of_segments.begin();
  std::list< std::pair<K::Segment_3, UInt> >::iterator end = list_of_segments.end();

  Array<Real>::vector_iterator nodes_it = nodes->begin(d);
  Array<UInt>::vector_iterator connectivity_it = connectivity.begin(2);

  for (; it != end ; ++it, nodes_it += 2, ++connectivity_it) {
    UInt node_id = nodes_it - nodes->begin(d);
    (*connectivity_it)(0) = node_id;
    (*connectivity_it)(1) = node_id + 1;

    for (UInt j = 0 ; j < d ; j++) {
      (*nodes_it)(j) = it->first.source()[j];
      (*(nodes_it + 1))(j) = it->first.target()[j];
    }
  }

  Mesh * mesh = new Mesh(d, *nodes, "interface_mesh");
  mesh->addConnectivityType(_segment_2);
  mesh->getConnectivity(_segment_2).copy(connectivity);

  if (associated_id) {
    delete associated_id;
    associated_id = NULL;
  }

  if (associated_type) {
    delete associated_type;
    associated_type = NULL;
  }

  associated_id = new Array<UInt>(mesh->getNbElement(_segment_2), 1, (UInt) 0);
  associated_type = new Array<Element>(mesh->getNbElement(_segment_2));

  Array<UInt>::scalar_iterator associated_id_it = associated_id->begin();
  Array<Element>::scalar_iterator associated_type_it = associated_type->begin();

  it = list_of_segments.begin();

  for (; it != end ; ++it, ++associated_id_it, ++associated_type_it) {
    *associated_id_it = it->second;
    associated_type_it->type = el_type;
  }

  mesh->registerData<UInt>("associated_id").setArray(_segment_2, _not_ghost, *associated_id);
  mesh->registerData<Element>("associated_type").setArray(_segment_2, _not_ghost, *associated_type);

  return mesh;
}

template<UInt d, ElementType el_type>
void MeshGeomFactory<d, el_type>::constructSegments(const std::list< typename TreeTypeHelper<d, el_type>::linear_intersection > & intersections,
                                                    std::list< std::pair<K::Segment_3, UInt> > & segments) {
  typename std::list<typename TreeTypeHelper<d, el_type>::linear_intersection>::const_iterator int_it = intersections.begin(),
                                                                                               int_end = intersections.end();

  for (; int_it != int_end ; ++int_it) {
    if (const K::Segment_3 * segment = boost::get<K::Segment_3>(&((*int_it)->first))) {
      std::pair<K::Segment_3, UInt> segment_id(*segment, (*int_it)->second);
      segments.push_back(segment_id);
    }
  }

  segments.unique(CompareSegmentPairs());
}



__END_AKANTU__

#endif // _AKANTU_MESH_GEOM_FACTORY_TMPL_HH__
