#===============================================================================
# @file   cgal.cmake
#
# @author Lucas Frerot <lucas.frerot@epfl.ch>
# @author Clement Roux <clement.roux@epfl.ch>
#
# @date creation: Thu Feb 19 2015
# @date last modification: Wed Jan 20 2016
#
# @brief  package description for CGAL
#
# @section LICENSE
#
# Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
# (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free  software: you can redistribute it and/or  modify it under the
# terms  of the  GNU Lesser  General Public  License as  published by  the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
# details.
#
# You should  have received  a copy  of the GNU  Lesser General  Public License
# along with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================
package_declare(CGAL EXTERNAL
  DESCRIPTION "Add CGAL support in akantu"
  )

package_is_activated(CGAL _is_activated)
package_on_enabled_script(CGAL
  "
  set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE
    CACHE BOOL \"Tells CGAL cmake to shut up\")
  set(CGAL_DISABLE_ROUNDING_MATH_CHECK ON
    CACHE BOOL \"Disable rounding math check in CGAL. This permits Valgrind to run.\")
  mark_as_advanced(
    CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE
    CGAL_DISABLE_ROUNDING_MATH_CHECK
    )
  ")

package_declare_sources(CGAL
  geometry/mesh_geom_common.hh

  geometry/mesh_geom_abstract.hh

  geometry/mesh_geom_factory.hh
  geometry/mesh_geom_factory_tmpl.hh

  geometry/mesh_abstract_intersector.hh
  geometry/mesh_abstract_intersector_tmpl.hh

  geometry/mesh_geom_intersector.hh
  geometry/mesh_geom_intersector_tmpl.hh

  geometry/mesh_segment_intersector.hh
  geometry/mesh_segment_intersector_tmpl.hh

  geometry/mesh_sphere_intersector.hh
  geometry/mesh_sphere_intersector_tmpl.hh

  geometry/tree_type_helper.hh
  geometry/geom_helper_functions.hh

  geometry/aabb_primitives/triangle.hh
  geometry/aabb_primitives/line_arc.hh
  geometry/aabb_primitives/tetrahedron.hh

  geometry/aabb_primitives/aabb_primitive.hh
  geometry/aabb_primitives/aabb_primitive.cc
  )

package_declare_documentation(CGAL
  "This package allows the use of CGAL's geometry algorithms in Akantu. Note that it needs a version of CGAL $\\geq$ 4.5 and needs activation of boost's system component."
  ""
  "CGAL checks with an assertion that the compilation flag \\shellcode{-frounding-math} is activated, which forbids the use of Valgrind on any code compilated with the package."
  )

package_set_package_system_dependency(CGAL deb-src "libcgal-dev >= 4.5")
