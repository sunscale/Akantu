/**
 * @file   mesh_geom_spherical_traits.hh
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Tue Jun 2 2015
 * @date last modification: Tue Jun 2 2015
 *
 * @brief Custom AABBGeomTraits to allow intersections in spherical kernel
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

#ifndef __AKANTU_MESH_GEOM_SPHERICAL_TRAITS_HH__
#define __AKANTU_MESH_GEOM_SPHERICAL_TRAITS_HH__

/* -------------------------------------------------------------------------- */

#include "mesh_geom_common.hh"

/* -------------------------------------------------------------------------- */

#define LINK_TYPEDEF(type) typedef typename SK::type type
#define TYPE_SEQ (Point_3)  \
                 (Sphere_3) \
                 (Construct_sphere_3) \
                 (Compute_closest_point_3)  \
                 (Has_on_bounded_side_3)    \
                 (Compute_squared_radius_3) \
                 (Compute_squared_distance_3)
#define AKANTU_LINK_TYPEDEF_MACRO(r, macro, elem) macro(elem);


__BEGIN_AKANTU__

/// Class used as AABBGeomTraits in the CGAL AABB tree algorithm
template<class SK>
class SphericalTraits {
  BOOST_PP_SEQ_FOR_EACH(AKANTU_LINK_TYPEDEF_MACRO, LINK_TYPEDEF, TYPE_SEQ);

};

__END_AKANTU__

#endif // __AKANTU_MESH_GEOM_SPHERICAL_KERNEL_HH__
