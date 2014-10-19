#===============================================================================
# @file   CMakeLists.txt
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
#
# @date creation: Fri Sep 03 2010
# @date last modification: Wed May 08 2013
#
# @brief  configuration for mesh partitioner tests
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

add_mesh(test_mesh_partitionate_mesh
  triangle.geo 2 2)

add_mesh(test_mesh_partitionate_mesh_data_mesh
  quad.geo 2 1)

register_test(test_mesh_partitionate_scotch
  SOURCES test_mesh_partitionate_scotch.cc
  DEPENDENCIES test_mesh_partitionate_mesh
  DIRECTORIES_TO_CREATE paraview
  )


register_test(test_mesh_partitionate_mesh_data
  SOURCES test_mesh_partitionate_mesh_data.cc
  DEPENDENCIES test_mesh_partitionate_mesh_data_mesh
  DIRECTORIES_TO_CREATE paraview)

