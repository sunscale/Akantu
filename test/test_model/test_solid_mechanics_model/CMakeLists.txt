#===============================================================================
# @file   CMakeLists.txt
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
#
# @date creation: Fri Sep 03 2010
# @date last modification: Tue Jan 30 2018
#
# @brief  configuratio for SolidMechanicsModel tests
#
# @section LICENSE
#
# Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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

add_akantu_test(test_materials "test_materials")
add_akantu_test(test_cohesive "cohesive_test")
add_akantu_test(test_embedded_interface "test_embedded_interface")
add_akantu_test(test_energies "test energies")
#===============================================================================

#===============================================================================
add_mesh(test_cube3d_two_mat_mesh cube_two_materials.geo 3 1)
register_test(test_solid_mechanics_model_reassign_material
  SOURCES test_solid_mechanics_model_reassign_material.cc
  DEPENDS test_cube3d_two_mat_mesh
  FILES_TO_COPY two_materials.dat
  PACKAGE parallel implicit
  PARALLEL
  )

#===============================================================================
register_test(test_solid_mechanics_model_material_eigenstrain
  SOURCES test_solid_mechanics_model_material_eigenstrain.cc
  FILES_TO_COPY cube_3d_tet_4.msh; material_elastic_plane_strain.dat
  PACKAGE implicit
  )

#===============================================================================
register_test(test_solid_mechanics_model_material_large_rotation
  SOURCES test_solid_mechanics_model_material_large_rotation.cc
  FILES_TO_COPY cube_3d_tet_4.msh; material_elastic.dat
  PACKAGE implicit
  )

#===============================================================================

register_test(test_material_selector
  SOURCES test_material_selector.cc
  FILES_TO_COPY material_selector.dat material_selector.msh
  PACKAGE core
  )

#===============================================================================
# dynamics tests
#===============================================================================
register_gtest_sources(
  SOURCES test_solid_mechanics_model_dynamics.cc
  FILES_TO_COPY test_solid_mechanics_model_dynamics_material.dat
  PACKAGE core
  )

register_gtest_test(test_solid_mechanics_model
  DEPENDS ${PATCH_TEST_BAR_MESHES}
  #bar_segment_2 bar_segment_3
  #bar_triangle_3 bar_triangle_6
  #bar_quadrangle_4 bar_quadrangle_8
  #bar_tetrahedron_4 bar_tetrahedron_10
  #bar_hexahedron_8 bar_hexahedron_20
  #bar_pentahedron_6 bar_pentahedron_15
  PARALLEL
  )
