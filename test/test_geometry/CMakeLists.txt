#===============================================================================
# @file   CMakeLists.txt
#
# @author Lucas Frerot <lucas.frerot@epfl.ch>
#
# @date creation: Fri Jun 19 2015
# @date last modification: Sat Jan 23 2016
#
# @brief  configuration for solver tests
#
# @section LICENSE
#
# Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free  software: you can redistribute it and/or  modify it under the terms  of the  GNU Lesser  General Public  License as published by  the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more details.
#
# You should  have received  a copy  of the GNU  Lesser General  Public License along with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
# @section DESCRIPTION
#
#===============================================================================

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
