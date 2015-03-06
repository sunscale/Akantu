/**
 * @file   mesh_geom_factory_tmpl.hh
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Thu Feb 26 2015
 * @date last modification: Fri Mar 6 2015
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
  primitive_list()
{}

template<UInt d, ElementType type>
MeshGeomFactory<d, type>::~MeshGeomFactory() {
  delete data_tree;
}

template<UInt d, ElementType type>
void MeshGeomFactory<d, type>::constructData() {
  AKANTU_DEBUG_IN();

  const GhostType ghost_type = _not_ghost;

  primitive_list.clear();

  UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);
  const Array<UInt> & connectivity = mesh.getConnectivity(type, ghost_type);
  const Array<Real> & nodes = mesh.getNodes();

  Array<UInt>::const_vector_iterator begin = connectivity.begin(nb_nodes_per_element);
  Array<UInt>::const_vector_iterator it    = connectivity.begin(nb_nodes_per_element);
  Array<UInt>::const_vector_iterator end   = connectivity.end(nb_nodes_per_element);

  /// This loop builds the list of primitives
  for (; it != end ; ++it) {
    const Vector<UInt> & el_connectivity = *it;
    Matrix<Real> node_coordinates(d, nb_nodes_per_element);

    for (UInt i = 0 ; i < nb_nodes_per_element ; i++)
      for (UInt j = 0 ; j < d ; j++)
        node_coordinates(j, i) = nodes(el_connectivity(i), j);

    addPrimitive(node_coordinates, it - begin);
  }

  delete data_tree;

  data_tree = new typename TreeTypeHelper<d, type>::tree(primitive_list.begin(), primitive_list.end());

  AKANTU_DEBUG_OUT();
}

// Need to implement this method for every case needed
template<>
void MeshGeomFactory<2, _triangle_3>::addPrimitive(const Matrix<Real> & node_coordinates, UInt id) {
  TreeTypeHelper<2, _triangle_3>::point_type a(node_coordinates(0, 0), node_coordinates(1, 0), 0.);
  TreeTypeHelper<2, _triangle_3>::point_type b(node_coordinates(0, 1), node_coordinates(1, 1), 0.);
  TreeTypeHelper<2, _triangle_3>::point_type c(node_coordinates(0, 2), node_coordinates(1, 2), 0.);

  Triangle<K> t(a, b, c);
  t.setId(id);
  primitive_list.push_back(t);
}

template<UInt d, ElementType el_type>
UInt MeshGeomFactory<d, el_type>::numberOfIntersectionsWithInterface(const K::Segment_3 & interface) const {
  return data_tree->number_of_intersected_primitives(interface);
}

template<UInt d, ElementType el_type>
Mesh * MeshGeomFactory<d, el_type>::meshOfLinearInterface(const std::pair<K::Segment_3, std::string> & pair) {
  AKANTU_DEBUG_IN();

  const K::Segment_3 & interface = pair.first;
  const std::string & id_str = pair.second;

  UInt number_of_intersections = this->numberOfIntersectionsWithInterface(interface);

  if (!number_of_intersections) {
    return NULL;
  }

  std::list<typename TreeTypeHelper<d, el_type>::linear_intersection> list_of_intersections;
  std::list< std::pair<K::Segment_3, UInt> > list_of_segments; // Contains no duplicate elements

  /// Compute all the intersection pairs (segment + element id) and remove duplicates
  data_tree->all_intersections(interface, std::back_inserter(list_of_intersections));
  this->constructSegments(list_of_intersections, list_of_segments);

  UInt nb_elements = list_of_segments.size();

  /// Temporary arrays for storing nodes and connectivity
  Array<Real> nodes(2 * nb_elements, d, "factory_nodes");
  Array<UInt> connectivity(nb_elements, 2, "factory_connectivity");

  /// Temporary arrays for storing associated element id and type
  Array<UInt> associated_id(nb_elements, 1, (UInt) 0, "factory_associated_id");
  Array<Element> associated_type(nb_elements, "factory_associated_type");

  std::list< std::pair<K::Segment_3, UInt> >::iterator it  = list_of_segments.begin();
  std::list< std::pair<K::Segment_3, UInt> >::iterator end = list_of_segments.end();

  Array<Real>::vector_iterator nodes_it = nodes.begin(d);
  Array<UInt>::vector_iterator connectivity_it = connectivity.begin(2);

  Array<UInt>::scalar_iterator associated_id_it = associated_id.begin();
  Array<Element>::scalar_iterator associated_type_it = associated_type.begin();

  /// Loop over the intersections pairs (segment, id)
  for (; it != end ; ++it,
                     nodes_it += 2,
                     ++connectivity_it,
                     ++associated_id_it,
                     ++associated_type_it) {
    UInt node_id = nodes_it - nodes.begin(d);

    /// Build element connectivity
    (*connectivity_it)(0) = node_id;
    (*connectivity_it)(1) = node_id + 1;

    /// Copy nodes
    for (UInt j = 0 ; j < d ; j++) {
      (*nodes_it)(j) = it->first.source()[j];
      (*(nodes_it + 1))(j) = it->first.target()[j];
    }

    /// Copy associated element info
    *associated_id_it = it->second;
    associated_type_it->type = el_type;
  }

  Mesh * interface_mesh = new Mesh(d, "factory_interface_mesh" + id_str);
  interface_mesh->getNodes().copy(nodes);

  interface_mesh->addConnectivityType(_segment_2);
  interface_mesh->getConnectivity(_segment_2).copy(connectivity);

  /// Register associated element info in interface_mesh
  interface_mesh->registerData<UInt>("associated_id").alloc(nb_elements, 1, _segment_2, _not_ghost);
  interface_mesh->registerData<Element>("associated_type").alloc(nb_elements, 1, _segment_2, _not_ghost);

  /// Copy associated element info
  interface_mesh->getData<UInt>("associated_id")(_segment_2, _not_ghost).copy(associated_id);
  interface_mesh->getData<Element>("associated_type")(_segment_2, _not_ghost).copy(associated_type);

  AKANTU_DEBUG_OUT();

  return interface_mesh;
}

template<UInt d, ElementType el_type>
void MeshGeomFactory<d, el_type>::constructSegments(const std::list< typename TreeTypeHelper<d, el_type>::linear_intersection > & intersections,
                                                    std::list< std::pair<K::Segment_3, UInt> > & segments) {
  AKANTU_DEBUG_IN();

  typename std::list<typename TreeTypeHelper<d, el_type>::linear_intersection>::const_iterator int_it = intersections.begin(),
                                                                                               int_end = intersections.end();

  for (; int_it != int_end ; ++int_it) {
    if (const K::Segment_3 * segment = boost::get<K::Segment_3>(&((*int_it)->first))) {
      std::pair<K::Segment_3, UInt> segment_id(*segment, (*int_it)->second);
      segments.push_back(segment_id);
    }
  }

  segments.unique(CompareSegmentPairs());

  AKANTU_DEBUG_OUT();
}



__END_AKANTU__

#endif // _AKANTU_MESH_GEOM_FACTORY_TMPL_HH__
