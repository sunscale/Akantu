/**
 * @file   tree_type_helper.hh
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Fri Feb 27 2015
 * @date last modification: Thu Mar 5 2015
 *
 * @brief  Converts element types of a mesh to CGAL primitive types
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

#ifndef __AKANTU_TREE_TYPE_HELPER_HH__
#define __AKANTU_TREE_TYPE_HELPER_HH__

#include "aka_common.hh"

#include "line_arc.hh"
#include "triangle.hh"
#include "tetrahedron.hh"

#include "aabb_primitive.hh"

#include "mesh_geom_common.hh"
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>

__BEGIN_AKANTU__
  
/* -------------------------------------------------------------------------- */

template<typename iterator>
struct VoidTree {
  VoidTree(const iterator & begin, const iterator & end) {}
};

/// Helper class used to ease the use of CGAL AABB tree algorithm
template<class Primitive, class Kernel>
struct TreeTypeHelper {
  static const bool is_valid = false;

  typedef Primitive primitive_type;
  typedef typename std::list<primitive_type> container_type;
  typedef typename container_type::iterator iterator;
  typedef typename container_type::const_iterator const_iterator;
  typedef typename CGAL::Point_3<Kernel> point_type;
  typedef VoidTree<iterator> tree;
};

/// Helper class used to ease the use of intersections
template<class TTHelper, class Query>
struct IntersectionTypeHelper;


/**
 * Macro used to specialize TreeTypeHelper
 * @param my_primitive associated primitive type
 * @param my_query query_type
 * @param my_kernel kernel type
 */
#define TREE_TYPE_HELPER_MACRO(my_primitive, my_query, my_kernel)      \
  template<>                                                            \
  struct TreeTypeHelper<my_primitive<my_kernel>, my_kernel> {            \
    static const bool is_valid = true;                                          \
    typedef my_primitive<my_kernel> primitive_type;                       \
    typedef my_primitive##_primitive aabb_primitive_type;                  \
    typedef CGAL::Point_3<my_kernel> point_type;                            \
                                                                             \
    typedef std::list<primitive_type> container_type;                         \
    typedef container_type::iterator iterator;                                 \
    typedef CGAL::AABB_traits<my_kernel, aabb_primitive_type> aabb_traits_type; \
    typedef CGAL::AABB_tree<aabb_traits_type> tree;                              \
    typedef tree::Primitive_id id_type;                                           \
  };                                                                               \
                                                                                    \
  template<>                                                                         \
  struct IntersectionTypeHelper<TreeTypeHelper<my_primitive<my_kernel>, my_kernel>, my_query> {             \
    typedef boost::optional<                                                                                 \
      TreeTypeHelper<my_primitive<my_kernel>, my_kernel>::tree::Intersection_and_primitive_id<my_query>::Type \
    > intersection_type;                                                                                       \
  }

TREE_TYPE_HELPER_MACRO(Triangle, Cartesian::Segment_3, Cartesian);
//TREE_TYPE_HELPER_MACRO(Line_arc, Spherical::Sphere_3, Spherical);

#undef TREE_TYPE_HELPER_MACRO

__END_AKANTU__

#endif // __AKANTU_TREE_TYPE_HELPER_HH__
