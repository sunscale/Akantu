#===============================================================================
# @file   CMakeLists.txt
#
# @author Lucas Frerot <lucas.frerot@epfl.ch>
#
# @date creation: Fri Feb 27 2015
# @date last modification: Fri Feb 27 2015
#
# @brief  configuration for solver tests
#
# @section LICENSE
#
# Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
# @section DESCRIPTION
#
#===============================================================================

add_mesh(test_geometry_intersection_triangle_3_mesh triangle_3.geo 2 1)
add_mesh(test_geometry_intersection_tetrahedron_4_mesh tetrahedron_4.geo 3 1)

register_test(test_geometry_predicates
  SOURCES test_geometry_predicates.cc
  PACKAGE CGAL
  )

register_test(test_geometry_intersection
  SOURCES test_geometry_intersection.cc
  FILES_TO_COPY test_geometry_triangle.msh
  PACKAGE CGAL
  )

register_test(test_segment_intersection_triangle_3
  SOURCES test_segment_intersection_triangle_3.cc
  FILES_TO_COPY test_geometry_triangle.msh
  PACKAGE CGAL
  )

register_test(test_segment_intersection_tetrahedron_4
  SOURCES test_segment_intersection_tetrahedron_4.cc
  FILES_TO_COPY test_geometry_tetrahedron.msh
  PACKAGE CGAL
  )
