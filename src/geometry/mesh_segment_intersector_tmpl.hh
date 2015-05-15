/**
 * @file mesh_segment_intersector_tmpl.hh
 *
 * @author Lucas Frerot
 *
 * @date creation: Wed Apr 29 2015
 * @date last modification: Wed Apr 29 2015
 *
 * @brief Computation of mesh intersection with segments
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2015 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef __AKANTU_MESH_SEGMENT_INTERSECTOR_TMPL_HH__
#define __AKANTU_MESH_SEGMENT_INTERSECTOR_TMPL_HH__

#include "aka_common.hh"
#include "tree_type_helper.hh"

#include <CGAL/Cartesian.h>

__BEGIN_AKANTU__

typedef CGAL::Cartesian<Real> K;

template<UInt dim, ElementType type>
MeshSegmentIntersector<dim, type>::MeshSegmentIntersector(const Mesh & mesh, Mesh & result_mesh):
  parent_type(mesh),
  result_mesh(result_mesh),
  current_physical_name()
{
  this->constructData();
}

template<UInt dim, ElementType type>
MeshSegmentIntersector<dim, type>::~MeshSegmentIntersector()
{}

template<UInt dim, ElementType type>
void MeshSegmentIntersector<dim, type>::computeIntersectionQuery(const K::Segment_3 & query) {
  AKANTU_DEBUG_IN();
  
  result_mesh.addConnectivityType(_segment_2, _not_ghost);
  result_mesh.addConnectivityType(_segment_2, _ghost);

  std::list<result_type> result_list;
  std::list<std::pair<K::Segment_3, UInt> > segment_list;

  this->factory.getTree().all_intersections(query, std::back_inserter(result_list));
  this->computeSegments(result_list, segment_list, query);

  // Arrays for storing nodes and connectivity
  Array<Real> & nodes = result_mesh.getNodes();
  Array<UInt> & connectivity = result_mesh.getConnectivity(_segment_2);

  // Arrays for storing associated element and physical name
  bool valid_elemental_data = true;
  Array<Element> * associated_element = NULL;
  Array<std::string> * associated_physical_name = NULL;

  try {
    associated_element = &result_mesh.getData<Element>("associated_element", _segment_2);
    associated_physical_name = &result_mesh.getData<std::string>("physical_names", _segment_2);
  } catch (debug::Exception & e) {
    valid_elemental_data = false;
  }

  std::list<pair_type>::iterator
    it = segment_list.begin(),
    end = segment_list.end();

  // Loop over the segment pairs
  for (; it != end ; ++it) {
    Vector<UInt> segment_connectivity(2);
    segment_connectivity(0) = result_mesh.getNbNodes();
    segment_connectivity(1) = result_mesh.getNbNodes() + 1;
    connectivity.push_back(segment_connectivity);

    // Copy nodes
    Vector<Real> source(dim), target(dim);
    for (UInt j = 0 ; j < dim ; j++) {
      source(j) = it->first.source()[j];
      target(j) = it->first.target()[j];
    }

    nodes.push_back(source);
    nodes.push_back(target);

    // Copy associated element info
    if (valid_elemental_data) {
      associated_element->push_back(Element(type, it->second));
      associated_physical_name->push_back(current_physical_name);
    }
  }

  AKANTU_DEBUG_OUT();
}

template<UInt dim, ElementType type>
void MeshSegmentIntersector<dim, type>::computeIntersectionQueryList(const std::list<K::Segment_3> & query_list,
                                                                     const std::string & physical_name) {
  AKANTU_DEBUG_IN();
  
  current_physical_name = physical_name;
  parent_type::computeIntersectionQueryList(query_list);
  
  AKANTU_DEBUG_OUT();
}

template<UInt dim, ElementType type>
void MeshSegmentIntersector<dim, type>::computeIntersectionQueryList(const std::list<K::Segment_3> & query_list) {
  AKANTU_DEBUG_IN();
  
  parent_type::computeIntersectionQueryList(query_list);
  
  AKANTU_DEBUG_OUT();
}

template<UInt dim, ElementType type>
void MeshSegmentIntersector<dim, type>::computeSegments(const std::list<result_type> & intersections,
                                                        std::list<pair_type> & segments,
                                                        const K::Segment_3 & query) {
  AKANTU_DEBUG_IN();
  
  /*
   * Number of intersections = 0 means
   *
   * - query is completely outside mesh
   * - query is completely inside primitive
   *
   * We try to determine the case and still construct the segment list
   */
  if (intersections.size() == 0) {
    // We look at all the primitives intersected by two rays
    // If there is one primitive in common, then query is inside
    // that primitive
    K::Ray_3 ray1(query.source(), query.target());
    K::Ray_3 ray2(query.target(), query.source());

    std::list<UInt> ray1_results, ray2_results;

    this->factory.getTree().all_intersected_primitives(ray1, std::back_inserter(ray1_results));
    this->factory.getTree().all_intersected_primitives(ray2, std::back_inserter(ray2_results));

    std::set<UInt> id_set;
    bool inside_primitive = false;
    UInt primitive_id = 0;

    std::list<UInt>::iterator
      ray1_it  = ray1_results.begin(),
      ray2_it  = ray2_results.begin(),
      ray1_end = ray1_results.end(),
      ray2_end = ray2_results.end();

    // Populate the set with elements of first list
    for (; ray1_it != ray1_end ; ++ray1_it)
      id_set.insert(id_set.begin(), *ray1_it);

    // Test if first list contains an element of second list
    for (; ray2_it != ray2_end && !inside_primitive ; ++ray2_it) {
      if (id_set.find(*ray2_it) != id_set.end()) {
        inside_primitive = true;
        primitive_id = *ray2_it;
      }
    }

    if (inside_primitive) {
      segments.push_back(std::make_pair(query, primitive_id));
    }
  }

  else {
    typename std::list<result_type>::const_iterator
      it = intersections.begin(),
         end = intersections.end();

    for(; it != end ; ++it) {
      UInt el = (*it)->second;

      // Result of intersection is a segment
      if (const K::Segment_3 * segment = boost::get<K::Segment_3>(&((*it)->first))) {
        // Check if the segment was alread created
        if (std::find_if(segments.begin(), segments.end(), IsSameSegment(*segment)) == segments.end()) {
          segments.push_back(std::make_pair(*segment, el));
        }
      }

      // Result of intersection is a point
      else if (const K::Point_3 * point = boost::get<K::Point_3>(&((*it)->first))) {
        // We only want to treat points differently if we're in 3D with Tetra4 elements
        // This should be optimized by compilator
        if (dim == 3 && type == _tetrahedron_4) {
          UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
          TreeTypeHelper<Triangle<K>, K>::container_type facets;

          const Array<Real> & nodes = this->mesh.getNodes();
          Array<UInt>::const_vector_iterator
            connectivity_vec = this->mesh.getConnectivity(type).begin(nb_nodes_per_element);

          const Vector<UInt> & el_connectivity = connectivity_vec[el];

          Matrix<Real> node_coordinates(dim, nb_nodes_per_element);
          for (UInt i = 0 ; i < nb_nodes_per_element ; i++)
            for (UInt j = 0 ; j < dim ; j++)
              node_coordinates(j, i) = nodes(el_connectivity(i), j);

          this->factory.addPrimitive(node_coordinates, el, facets);

          // Local tree
          TreeTypeHelper<Triangle<K>, K>::tree * local_tree =
            new TreeTypeHelper<Triangle<K>, K>::tree(facets.begin(), facets.end());

          // Compute local intersections (with current element)
          std::list<result_type> local_intersections;
          local_tree->all_intersections(query, std::back_inserter(local_intersections));

          if (local_intersections.size() == 1) {
            TreeTypeHelper<Triangle<K>, K>::point_type
              a(node_coordinates(0, 0), node_coordinates(1, 0), node_coordinates(2, 0)),
              b(node_coordinates(0, 1), node_coordinates(1, 1), node_coordinates(2, 1)),
              c(node_coordinates(0, 2), node_coordinates(1, 2), node_coordinates(2, 2)),
              d(node_coordinates(0, 3), node_coordinates(1, 3), node_coordinates(2, 3));
            K::Tetrahedron_3 tetra(a, b, c, d);
            K::Point_3 inside_point =
              (tetra.has_on_bounded_side(query.source())) ? query.source() : query.target();
            K::Segment_3 seg(inside_point, *point);
            segments.push_back(std::make_pair(seg, el));
          }

          else {
            typename std::list<result_type>::const_iterator
              local_it = local_intersections.begin(),
                       local_end = local_intersections.end();

            for (; local_it != local_end ; ++local_it) {
              if (const K::Point_3 * local_point = boost::get<K::Point_3>(&((*local_it)->first))) {
                if (!comparePoints(*point, *local_point)) {
                  K::Segment_3 seg(*point, *local_point);
                  if(std::find_if(segments.begin(), segments.end(), IsSameSegment(seg)) == segments.end())
                    segments.push_back(std::make_pair(seg, el));
                }
              }
            }
          }

          delete local_tree;
        }
      }
    }
  }
  AKANTU_DEBUG_OUT();
}

__END_AKANTU__

#endif // __AKANTU_MESH_SEGMENT_INTERSECTOR_TMPL_HH__
