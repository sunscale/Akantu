#  @file   CMakeLists.txt
#  @author Fabian Barras <fabian.barras@epfl.ch>
#  @date   Wed Jan 22 15:09:31 2014
# 
#  @section LICENSE
# 
#  Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
#  Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
# 
#  Akantu is free  software: you can redistribute it and/or  modify it under the
#  terms  of the  GNU Lesser  General Public  License as  published by  the Free
#  Software Foundation, either version 3 of the License, or (at your option) any
#  later version.
# 
#  Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
#  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
#  A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
#  details.
# 
#  You should  have received  a copy  of the GNU  Lesser General  Public License
#  along with Akantu. If not, see <http://www.gnu.org/licenses/>.
# 
#===============================================================================

#===============================================================================
add_mesh(mesh_dcb_2d mesh_dcb_2d.geo 2 2)

register_example(test_dcb_2d
  test_dcb_2d.cc)

add_dependencies(test_dcb_2d
  mesh_dcb_2d)

#===============================================================================
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/input_test.dat
  ${CMAKE_CURRENT_BINARY_DIR}/input_test.dat
  COPYONLY
  )

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview)