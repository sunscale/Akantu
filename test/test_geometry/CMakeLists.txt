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

add_mesh(test_geometry_intersection_mesh mesh.geo 2 1)
register_test(test_geometry_intersection
  SOURCES test_geometry_intersection.cc
  DEPENDENCIES test_geometry_intersection_mesh
  PACKAGE CGAL
  )
