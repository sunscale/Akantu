#===============================================================================
# @file   CMakeLists.txt
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Sep 03 2010
# @date last modification: Fri May 03 2013
#
# @brief  configuration for MeshIO tests
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

add_mesh(test_msh_cube cube.geo 3 1)
add_mesh(test_msh_cube_physical_names cube_physical_names.geo 3 1)

register_test(test_mesh_io_msh
  SOURCES test_mesh_io_msh.cc
  DEPENDENCIES test_msh_cube
  DIRECTORIES_TO_CREATE paraview
  )
  
register_test(test_mesh_io_msh_physical_names
  SOURCES test_mesh_io_msh_physical_names.cc
  DEPENDENCIES test_msh_cube_physical_names
)

#===============================================================================
#register_test(test_mesh_io_diana test_mesh_io_diana.cc)
#copy_files(test_mesh_io_diana dam.dat)
#file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview/mesh_io_diana)

