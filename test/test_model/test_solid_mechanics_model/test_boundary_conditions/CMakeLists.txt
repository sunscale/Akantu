#===============================================================================
# @file   CMakeLists.txt
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
#
# @date creation: Fri Nov 26 2010
# @date last modification: Fri Feb 14 2014
#
# @brief  configuration for materials tests
#
# @section LICENSE
#
# Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

add_mesh(periodic_plate periodic_plate.geo 2 1)
register_test(test_node_selection
  SOURCES  test_node_selection.cc 
  DEPENDS periodic_plate
  FILES_TO_COPY material.dat
  DIRECTORIES_TO_CREATE paraview
  PACKAGE implicit extra_materials
  )
