/**
 * @file   tree_type_helper.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Fri Jan 04 2013
 * @date last modification: Thu Feb 01 2018
 *
 * @brief  Converts element types of a mesh to CGAL primitive types
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
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
#include "tetrahedron.hh"
#include "triangle.hh"

#include "aabb_primitive.hh"

#include "mesh_geom_common.hh"
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>

namespace akantu {

/* -------------------------------------------------------------------------- */

/// Replacement class for algorithm that can't use the AABB tree types
template <typename iterator> struct VoidTree {
  VoidTree(const iterator & /*begin*/, const iterator & /*end*/) {}
};

/// Helper class used to ease the use of CGAL AABB tree algorithm
template <class Primitive, class Kernel> struct TreeTypeHelper {
  static const bool is_valid = false;

  typedef Primitive primitive_type;
  typedef typename std::list<primitive_type> container_type;
  typedef typename container_type::iterator iterator;
  typedef typename container_type::const_iterator const_iterator;
  typedef typename CGAL::Point_3<Kernel> point_type;
  typedef VoidTree<iterator> tree;
};

/// Helper class used to ease the use of intersections
template <class TTHelper, class Query> struct IntersectionTypeHelper;

/**
 * Macro used to specialize TreeTypeHelper
 * @param my_primitive associated primitive type
 * @param my_query query_type
 * @param my_kernel kernel type
 */
#define TREE_TYPE_HELPER_MACRO(my_primitive, my_query, my_kernel)              \
  template <> struct TreeTypeHelper<my_primitive<my_kernel>, my_kernel> {      \
    static const bool is_valid = true;                                         \
    typedef my_primitive<my_kernel> primitive_type;                            \
    typedef my_primitive##_primitive aabb_primitive_type;                      \
    typedef CGAL::Point_3<my_kernel> point_type;                               \
                                                                               \
    typedef std::list<primitive_type> container_type;                          \
    typedef container_type::iterator iterator;                                 \
    typedef CGAL::AABB_traits<my_kernel, aabb_primitive_type>                  \
        aabb_traits_type;                                                      \
    typedef CGAL::AABB_tree<aabb_traits_type> tree;                            \
    typedef tree::Primitive_id id_type;                                        \
  };                                                                           \
                                                                               \
  template <>                                                                  \
  struct IntersectionTypeHelper<                                               \
      TreeTypeHelper<my_primitive<my_kernel>, my_kernel>, my_query> {          \
    typedef boost::optional<TreeTypeHelper<                                    \
        my_primitive<my_kernel>,                                               \
        my_kernel>::tree::Intersection_and_primitive_id<my_query>::Type>       \
        intersection_type;                                                     \
  }

TREE_TYPE_HELPER_MACRO(Triangle, cgal::Cartesian::Segment_3, cgal::Cartesian);
// TREE_TYPE_HELPER_MACRO(Line_arc, cgal::Spherical::Sphere_3, cgal::Spherical);

#undef TREE_TYPE_HELPER_MACRO

} // namespace akantu

#endif // __AKANTU_TREE_TYPE_HELPER_HH__
