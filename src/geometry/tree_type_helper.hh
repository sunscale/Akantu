/**
 * @file   tree_type_helper.hh
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Fri Feb 27 2015
 * @date last modification: Fri Feb 27 2015
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

#include <CGAL/Cartesian.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>

__BEGIN_AKANTU__
  
typedef CGAL::Cartesian<Real> Kernel;

template<UInt d, ElementType el_type>
struct TreeTypeHelper;

template<>
struct TreeTypeHelper<2, _triangle_3> {
    typedef Kernel::Triangle_3 primitive_type;
    typedef Kernel::Point_3 point_type;

    typedef std::list<primitive_type> container_type;
    typedef container_type::iterator iterator;
    typedef CGAL::AABB_triangle_primitive<Kernel, iterator> aabb_primitive_type;
    typedef CGAL::AABB_traits<Kernel, aabb_primitive_type> aabb_traits_type;
    typedef CGAL::AABB_tree<aabb_traits_type> tree;
};

__END_AKANTU__

#endif // __AKANTU_TREE_TYPE_HELPER_HH__
