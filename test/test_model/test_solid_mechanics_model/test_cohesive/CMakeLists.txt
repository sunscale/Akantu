#===============================================================================
# @file   CMakeLists.txt
#
# @author Marco Vocialta <marco.vocialta@epfl.ch>
#
# @date creation: Fri Oct 22 2010
# @date last modification: Wed Feb 21 2018
#
# @brief  configuration for cohesive elements tests
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

#add_akantu_test(test_cohesive_intrinsic "test_cohesive_intrinsic")
#add_akantu_test(test_cohesive_extrinsic "test_cohesive_extrinsic")
#add_akantu_test(test_cohesive_intrinsic_impl "test_cohesive_intrinsic_impl")
#add_akantu_test(test_cohesive_1d_element "test_cohesive_1d_element")
#add_akantu_test(test_cohesive_extrinsic_implicit "test_cohesive_extrinsic_implicit")

add_akantu_test(test_materials "test_cohesive_materials")

add_akantu_test(test_cohesive_buildfragments "test_cohesive_buildfragments")
add_akantu_test(test_cohesive_insertion "test_cohesive_insertion")
add_akantu_test(test_cohesive_linear_friction "test_cohesive_linear_friction")
#add_akantu_test(test_parallel_cohesive "parallel_cohesive_test")

set(_meshes)

add_mesh(cohesive_1d_2_seg data/cohesive_1D.geo
  DIM 2 ORDER 1
  OUTPUT _cohesive_1d_2_segment_2.msh)
list(APPEND _meshes cohesive_1d_2_seg)

add_mesh(cohesive_2d_4_tri data/cohesive_strait_2D.geo
  DIM 2 ORDER 1
  OUTPUT _cohesive_2d_4_triangle_3.msh)
list(APPEND _meshes cohesive_2d_4_tri)

add_mesh(cohesive_2d_6_tri data/cohesive_strait_2D.geo
  DIM 2 ORDER 2
  OUTPUT _cohesive_2d_6_triangle_6.msh)
list(APPEND _meshes cohesive_2d_6_tri)

add_mesh(cohesive_2d_4_quad data/cohesive_strait_2D_structured.geo
  DIM 2 ORDER 1
  OUTPUT _cohesive_2d_4_quadrangle_4.msh)
list(APPEND _meshes cohesive_2d_4_quad)

add_mesh(cohesive_2d_6_quad data/cohesive_strait_2D_structured.geo
  DIM 2 ORDER 2
  OUTPUT _cohesive_2d_6_quadrangle_8.msh)
list(APPEND _meshes cohesive_2d_6_quad)

add_mesh(cohesive_2d_4_tri_quad data/cohesive_strait_2D_mixte.geo
  DIM 2 ORDER 1
  OUTPUT _cohesive_2d_4_triangle_3_quadrangle_4.msh)
list(APPEND _meshes cohesive_2d_4_tri_quad)

add_mesh(cohesive_2d_6_tri_quad data/cohesive_strait_2D_mixte.geo
  DIM 2 ORDER 2
  OUTPUT _cohesive_2d_6_triangle_6_quadrangle_8.msh)
list(APPEND _meshes cohesive_2d_6_tri_quad)

add_mesh(cohesive_3d_6_tet data/cohesive_strait_3D.geo
  DIM 3 ORDER 1
  OUTPUT _cohesive_3d_6_tetrahedron_4.msh)
list(APPEND _meshes cohesive_3d_6_tet)

add_mesh(cohesive_3d_12_tet data/cohesive_strait_3D.geo
  DIM 3 ORDER 2
  OUTPUT _cohesive_3d_12_tetrahedron_10.msh)
list(APPEND _meshes cohesive_3d_12_tet)

add_mesh(cohesive_3d_8_hex data/cohesive_strait_3D_structured.geo
  DIM 3 ORDER 1
  OUTPUT _cohesive_3d_8_hexahedron_8.msh)
list(APPEND _meshes cohesive_3d_8_hex)

add_mesh(cohesive_3d_16_hex data/cohesive_strait_3D_structured.geo
  DIM 3 ORDER 2
  OUTPUT _cohesive_3d_16_hexahedron_20.msh)
list(APPEND _meshes cohesive_3d_16_hex)

register_gtest_sources(
  SOURCES test_cohesive.cc
  PACKAGE cohesive_element
  DEPENDS ${_meshes}
  FILES_TO_COPY material_0.dat material_1.dat material_0_finite_def.dat
  )

register_gtest_test(test_solid_mechanics_model_cohesive
  PARALLEL
  PARALLEL_LEVEL 1 2
  )
# ==============================================================================
