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

#include <algorithm>

#include <CGAL/Cartesian.h>

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

typedef CGAL::Cartesian<Real> K;

template<UInt dim, ElementType type>
MeshGeomFactory<dim, type>::MeshGeomFactory(const Mesh & mesh) :
  MeshGeomAbstract(mesh),
  data_tree(NULL),
  primitive_list()
{}

template<UInt dim, ElementType type>
MeshGeomFactory<dim, type>::~MeshGeomFactory() {
  delete data_tree;
}

template<UInt dim, ElementType type>
void MeshGeomFactory<dim, type>::constructData() {
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
    Matrix<Real> node_coordinates(dim, nb_nodes_per_element);

    for (UInt i = 0 ; i < nb_nodes_per_element ; i++)
      for (UInt j = 0 ; j < dim ; j++)
        node_coordinates(j, i) = nodes(el_connectivity(i), j);

    addPrimitive(node_coordinates, it - begin);
  }

  delete data_tree;

  data_tree = new typename TreeTypeHelper<dim, type>::tree(primitive_list.begin(), primitive_list.end());

  AKANTU_DEBUG_OUT();
}

// 2D and _triangle_3 implementation
template<>
void MeshGeomFactory<2, _triangle_3>::addPrimitive(const Matrix<Real> & node_coordinates, UInt id) {
  TreeTypeHelper<2, _triangle_3>::point_type a(node_coordinates(0, 0), node_coordinates(1, 0), 0.);
  TreeTypeHelper<2, _triangle_3>::point_type b(node_coordinates(0, 1), node_coordinates(1, 1), 0.);
  TreeTypeHelper<2, _triangle_3>::point_type c(node_coordinates(0, 2), node_coordinates(1, 2), 0.);

  Triangle<K> t(a, b, c);
  t.setId(id);
  primitive_list.push_back(t);
}

// 3D and _tetrahedron_4 implementation
template<>
void MeshGeomFactory<3, _tetrahedron_4>::addPrimitive(const Matrix<Real> & node_coordinates, UInt id) {
  TreeTypeHelper<3, _tetrahedron_4>::point_type
    a(node_coordinates(0, 0), node_coordinates(1, 0), node_coordinates(2, 0)),
    b(node_coordinates(0, 1), node_coordinates(1, 1), node_coordinates(2, 1)),
    c(node_coordinates(0, 2), node_coordinates(1, 2), node_coordinates(2, 2)),
    d(node_coordinates(0, 3), node_coordinates(1, 3), node_coordinates(2, 3));

  Triangle<K>
    t1(a, b, c),
    t2(b, c, d),
    t3(c, d, a),
    t4(d, a, b);

  t1.setId(id);
  t2.setId(id);
  t3.setId(id);
  t4.setId(id);

  primitive_list.push_back(t1);
  primitive_list.push_back(t2);
  primitive_list.push_back(t3);
  primitive_list.push_back(t4);
}

template<UInt dim, ElementType el_type>
UInt MeshGeomFactory<dim, el_type>::numberOfIntersectionsWithInterface(const K::Segment_3 & interface) const {
  return data_tree->number_of_intersected_primitives(interface);
}

template<UInt dim, ElementType el_type>
void MeshGeomFactory<dim, el_type>::meshOfLinearInterface(const Interface & interface_pair, Mesh & interface_mesh) {
  AKANTU_DEBUG_IN();

  const K::Segment_3 & interface = interface_pair.first;

  UInt number_of_intersections = this->numberOfIntersectionsWithInterface(interface);

  if (!number_of_intersections) {
    AKANTU_DEBUG_WARNING("No intersection with interface");
    return;
  }

  std::list<typename TreeTypeHelper<dim, el_type>::linear_intersection> list_of_intersections;
  /// TODO change this list to a set using lessSegmentPair for optimization
  std::list< std::pair<K::Segment_3, UInt> > list_of_segments; // Contains no duplicate elements

  /// Compute all the intersection pairs (segment + element id) and remove duplicates
  data_tree->all_intersections(interface, std::back_inserter(list_of_intersections));
  this->constructSegments(list_of_intersections, list_of_segments, interface);

  /// Arrays for storing nodes and connectivity
  Array<Real> & nodes = interface_mesh.getNodes();
  Array<UInt> & connectivity = interface_mesh.getConnectivity(_segment_2);

  /// Arrays for storing associated element id and type
  Array<Element> & associated_element = interface_mesh.getData<Element>("associated_element", _segment_2);
  Array<std::string> & associated_material = interface_mesh.getData<std::string>("material", _segment_2);

  std::list<std::pair<K::Segment_3, UInt> >::iterator it  = list_of_segments.begin();
  std::list<std::pair<K::Segment_3, UInt> >::iterator end = list_of_segments.end();

  /// Loop over the intersections pairs (segment, id)
  for (; it != end ; ++it) {
    Vector<UInt> segment_connectivity(2);
    segment_connectivity(0) = interface_mesh.getNbNodes();
    segment_connectivity(1) = interface_mesh.getNbNodes() + 1;
    connectivity.push_back(segment_connectivity);

    /// Copy nodes
    Vector<Real> source(dim), target(dim);
    for (UInt j = 0 ; j < dim ; j++) {
      source(j) = it->first.source()[j];
      target(j) = it->first.target()[j];
    }

    nodes.push_back(source);
    nodes.push_back(target);

    /// Copy associated element info
    associated_element.push_back(Element(el_type, it->second));
    associated_material.push_back(interface_pair.second);
  }

  AKANTU_DEBUG_OUT();
}

template<UInt dim, ElementType el_type>
void MeshGeomFactory<dim, el_type>::constructSegments(
    const std::list< typename TreeTypeHelper<dim, el_type>::linear_intersection > & intersections,
    std::list<std::pair<K::Segment_3, UInt> > & segments,
    const K::Segment_3 & interface)
{
  AKANTU_DEBUG_IN();

  typename std::list<typename TreeTypeHelper<dim, el_type>::linear_intersection>::const_iterator
    int_it = intersections.begin(),
    int_end = intersections.end();

  for (; int_it != int_end ; ++int_it) {
    UInt el = (*int_it)->second;

    if (const K::Segment_3 * segment = boost::get<K::Segment_3>(&((*int_it)->first))) {
      if (std::find_if(segments.begin(), segments.end(), IsSameSegment(*segment)) == segments.end()) {
        // Use min() and max() for ordering
        K::Segment_3 seg_ord(segment->min(), segment->max());
        segments.push_back(std::make_pair(seg_ord, (*int_it)->second));
      }
    }

    else if (const K::Point_3 * point = boost::get<K::Point_3>(&((*int_it)->first))) {
      // We only want to treat points differently if we're in 3D with Tetra4 elements
      // This should be optimized by compilator
      if (dim == 3 && el_type == _tetrahedron_4) {
        UInt nb_facets = Mesh::getNbFacetsPerElement(el_type);
        typename TreeTypeHelper<dim, el_type>::container_type facets(nb_facets);

        // TODO Use mesh facets instead of this (should make O(n²) go to O(n))
        std::remove_copy_if(primitive_list.begin(),
                            primitive_list.end(),
                            facets.begin(),
                            BelongsNotToElement<dim, el_type>(el));

        typename TreeTypeHelper<dim, el_type>::tree * local_tree = 
          new typename TreeTypeHelper<dim, el_type>::tree(facets.begin(), facets.end());

        std::list< typename TreeTypeHelper<dim, el_type>::linear_intersection>
          local_intersections;

        local_tree->all_intersections(interface, std::back_inserter(local_intersections));

        typename std::list<typename TreeTypeHelper<dim, el_type>::linear_intersection>::const_iterator
          local_it = local_intersections.begin(),
          local_end = local_intersections.end();

        for (; local_it != local_end ; ++local_it) {
          if (const K::Point_3 * local_point = boost::get<K::Point_3>(&((*local_it)->first))) {
            if (!comparePoints(*point, *local_point)) {
              K::Segment_3 seg(*point, *local_point);
              K::Segment_3 seg_ord(seg.min(), seg.max());
              if (std::find_if(segments.begin(), segments.end(), IsSameSegment(seg_ord)) == segments.end())
                segments.push_back(std::make_pair(seg_ord, el));
            }
          }
        }

        delete local_tree;
      }
    }
  }

  AKANTU_DEBUG_OUT();
}



__END_AKANTU__

#endif // _AKANTU_MESH_GEOM_FACTORY_TMPL_HH__
