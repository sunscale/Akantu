#===============================================================================
# @file   CMakeLists.txt
#
# @author Fabian Barras <fabian.barras@epfl.ch>
#
# @date creation: Mon Dec 14 09:07:44 2015
# 
# @brief  Tests insertion of cohesive elements
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
add_mesh(swiss_cheese swiss_cheese.geo 2 2)

register_example(example_parser
  example_parser.cc
  USE IOHelper
  DEPENDS swiss_cheese
  DIRECTORIES_TO_CREATE paraview
  )

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/input_file.dat
  ${CMAKE_CURRENT_BINARY_DIR}/input_file.dat
  COPYONLY
  )
