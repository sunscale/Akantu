#===============================================================================
# @file   CMakeLists.txt
#
# @author Marco Vocialta <marco.vocialta@epfl.ch>
#
# @date creation: Wed Nov 05 2014
# @date last modification: Fri Oct 13 2017
#
# @brief  configuration for parallel test for extrinsic cohesive elements
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

add_mesh(test_cohesive_parallel_extrinsic_mesh mesh.geo 2 2)
add_mesh(test_cohesive_parallel_extrinsic_tetrahedron_mesh tetrahedron.geo 3 2)

register_test(test_cohesive_parallel_extrinsic
  SOURCES test_cohesive_parallel_extrinsic.cc
  DEPENDS test_cohesive_parallel_extrinsic_mesh
  PACKAGE parallel_cohesive_element
  FILES_TO_COPY material.dat
  DIRECTORIES_TO_CREATE paraview)

register_test(test_cohesive_parallel_extrinsic_tetrahedron
  SOURCES test_cohesive_parallel_extrinsic_tetrahedron.cc
  DEPENDS test_cohesive_parallel_extrinsic_tetrahedron_mesh
  PACKAGE parallel_cohesive_element
  FILES_TO_COPY material.dat
  DIRECTORIES_TO_CREATE paraview)

register_test(test_cohesive_parallel_extrinsic_tetrahedron_displacement
  SOURCES test_cohesive_parallel_extrinsic_tetrahedron_displacement.cc
  DEPENDS test_cohesive_parallel_extrinsic_tetrahedron_mesh
  PACKAGE parallel_cohesive_element
  FILES_TO_COPY material.dat
  DIRECTORIES_TO_CREATE paraview displacement)
