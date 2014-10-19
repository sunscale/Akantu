#===============================================================================
# @file   CMakeLists.txt
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
#
# @date creation: Sun May 01 2011
# @date last modification: Tue Nov 06 2012
#
# @brief  configuration for heat transfer model tests
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

#===============================================================================
add_mesh(test_heat_cube3d_mesh1 cube.geo 3 1 OUTPUT cube1.msh)
add_mesh(test_heat_square_mesh1 square.geo 2 1 OUTPUT square1.msh)

add_mesh(test_heat_line_mesh line.geo 1 1 OUTPUT line.msh)
register_test(test_heat_transfer_model_cube3d
  SOURCES test_heat_transfer_model_cube3d.cc
  DEPENDENCIES test_heat_cube3d_mesh1 test_heat_square_mesh1 test_heat_line_mesh
  FILES_TO_COPY material.dat
  DIRECTORIES_TO_CREATE paraview
  )

add_mesh(test_heat_cube_tet4 cube_tet.geo 3 1 OUTPUT cube_tet4.msh)
register_test(test_heat_transfer_model_cube3d_pbc
  SOURCES test_heat_transfer_model_cube3d_pbc.cc
  DEPENDENCIES test_heat_cube_tet4
  FILES_TO_COPY material.dat
  DIRECTORIES_TO_CREATE paraview
  )

add_mesh(test_heat_square_tri3 square_tri.geo 2 1 OUTPUT square_tri3.msh)
register_test(test_heat_transfer_model_square2d_pbc
  SOURCES test_heat_transfer_model_square2d_pbc.cc
  DEPENDENCIES test_heat_square_tri3
  FILES_TO_COPY material.dat
  DIRECTORIES_TO_CREATE paraview
  )

register_test(test_heat_transfer_model_square2d
  SOURCES test_heat_transfer_model_square2d.cc
  DEPENDENCIES test_heat_square_tri3
  FILES_TO_COPY material.dat
  DIRECTORIES_TO_CREATE paraview
  )


register_test(test_heat_transfer_model_cube3d_istropic_conductivity
  SOURCES test_heat_transfer_model_cube3d_istropic_conductivity.cc
  DEPENDENCIES test_heat_cube3d_mesh1
  FILES_TO_COPY material.dat
  DIRECTORIES_TO_CREATE paraview
  )

