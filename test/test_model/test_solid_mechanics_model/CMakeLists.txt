#===============================================================================
# @file   CMakeLists.txt
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
#
# @date creation: Fri Sep 03 2010
# @date last modification: Tue Jan 19 2016
#
# @brief  configuratio for SolidMechanicsModel tests
#
# @section LICENSE
#
# Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
# Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
# Solides)
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

add_akantu_test(test_materials "test_materials")
add_akantu_test(patch_tests "patch_tests")
add_akantu_test(test_cohesive "cohesive_test")
add_akantu_test(test_embedded_interface "test_embedded_interface")

#===============================================================================
add_mesh(test_solid_mechanics_model_square_mesh square.geo 2 1)
add_mesh(test_solid_mechanics_model_circle_mesh1 circle.geo 2 1 OUTPUT circle1.msh)
add_mesh(test_solid_mechanics_model_circle_mesh2 circle.geo 2 2 OUTPUT circle2.msh)

register_test(test_solid_mechanics_model_square
  SOURCES test_solid_mechanics_model_square.cc
  DEPENDS test_solid_mechanics_model_square_mesh
  FILES_TO_COPY material.dat test_cst_energy.pl
  DIRECTORIES_TO_CREATE paraview
  PACKAGE core
  )

register_test(test_solid_mechanics_model_circle_2
  SOURCES test_solid_mechanics_model_circle_2.cc
  DEPENDS test_solid_mechanics_model_circle_mesh2
  FILES_TO_COPY material.dat
  DIRECTORIES_TO_CREATE paraview
  PACKAGE core
  )

#===============================================================================
add_mesh(test_bar_traction_2d_mesh1 bar.geo 2 1 OUTPUT bar1.msh)
add_mesh(test_bar_traction_2d_mesh2 bar.geo 2 2 OUTPUT bar2.msh)
add_mesh(test_bar_traction_2d_mesh_structured1 bar_structured.geo 2 1 OUTPUT bar_structured1.msh)

register_test(test_solid_mechanics_model_bar_traction2d
  SOURCES test_solid_mechanics_model_bar_traction2d.cc
  DEPENDS test_bar_traction_2d_mesh1 test_bar_traction_2d_mesh2
  FILES_TO_COPY material.dat test_cst_energy.pl
  DIRECTORIES_TO_CREATE paraview
  PACKAGE core
  )

register_test(test_solid_mechanics_model_bar_traction2d_structured
  SOURCES test_solid_mechanics_model_bar_traction2d_structured.cc
  DEPENDS test_bar_traction_2d_mesh_structured1
  FILES_TO_COPY material.dat test_cst_energy.pl
  DIRECTORIES_TO_CREATE paraview
  PACKAGE core
  )

#===============================================================================
add_mesh(test_solid_mechanics_model_segment_mesh segment.geo 1 2)

register_test(test_solid_mechanics_model_bar_traction2d_parallel
  SOURCES test_solid_mechanics_model_bar_traction2d_parallel.cc
  DEPENDS test_bar_traction_2d_mesh2
  FILES_TO_COPY material.dat test_cst_energy.pl
  DIRECTORIES_TO_CREATE paraview
  PACKAGE parallel
  )

register_test(test_solid_mechanics_model_segment_parallel
  SOURCES test_solid_mechanics_model_segment_parallel.cc
  DEPENDS test_solid_mechanics_model_segment_mesh
  FILES_TO_COPY material.dat test_cst_energy.pl
  DIRECTORIES_TO_CREATE paraview
  PACKAGE parallel
  )

#===============================================================================
#register_test(test_solid_mechanics_model_bar_traction2d_mass_not_lumped
#  SOURCES test_solid_mechanics_model_bar_traction2d_mass_not_lumped.cc
#  DEPENDS test_bar_traction_2d_mesh1 test_bar_traction_2d_mesh2
#  FILES_TO_COPY material.dat
#  DIRECTORIES_TO_CREATE paraview
#  PACKAGE implicit
#  )

#===============================================================================
add_mesh(test_solid_mechanics_model_segment_mesh1 segment.geo 1 1 OUTPUT segment1.msh)
add_mesh(test_implicit_mesh1 square_implicit.geo 2 1 OUTPUT square_implicit1.msh)
add_mesh(test_implicit_mesh2 square_implicit.geo 2 2 OUTPUT square_implicit2.msh)

register_test(test_solid_mechanics_model_implicit_1d
  SOURCES test_solid_mechanics_model_implicit_1d.cc
  DEPENDS test_solid_mechanics_model_segment_mesh1
  FILES_TO_COPY material.dat
  DIRECTORIES_TO_CREATE paraview
  PACKAGE implicit
  )

register_test(test_solid_mechanics_model_implicit_2d
  SOURCES test_solid_mechanics_model_implicit_2d.cc
  DEPENDS test_implicit_mesh1 test_implicit_mesh2
  FILES_TO_COPY material.dat
  DIRECTORIES_TO_CREATE paraview
  PACKAGE implicit
  )

#===============================================================================
add_mesh(test_implicit_beam_2d_1 beam_2d.geo 2 1 OUTPUT beam_2d_lin.msh)
add_mesh(test_implicit_beam_2d_2 beam_2d.geo 2 2 OUTPUT beam_2d_quad.msh)
add_mesh(test_implicit_beam_3d_2 beam_3d.geo 3 2 OUTPUT beam_3d_quad.msh)
add_mesh(test_implicit_beam_3d_1 beam_3d.geo 3 1 OUTPUT beam_3d_lin.msh)
register_test(test_solid_mechanics_model_implicit_dynamic_2d
  SOURCES test_solid_mechanics_model_implicit_dynamic_2d.cc
  DEPENDS test_implicit_beam_2d_1 test_implicit_beam_2d_2 test_implicit_beam_3d_2 test_implicit_beam_3d_1
  FILES_TO_COPY  material_implicit_dynamic.dat
  DIRECTORIES_TO_CREATE paraview
  PACKAGE implicit
  )

#===============================================================================
# register_test(test_solid_mechanics_model_bar_traction2d_structured_pbc
#   SOURCES test_solid_mechanics_model_bar_traction2d_structured_pbc.cc
#   DEPENDS test_bar_traction_2d_mesh_structured1
#   FILES_TO_COPY material.dat test_cst_energy.pl
#   DIRECTORIES_TO_CREATE paraview
#   PACKAGE core
#   )

#add_mesh(test_pbc_parallel_mesh square_structured.geo 2 1 OUTPUT square_structured.msh)
#register_test(test_solid_mechanics_model_pbc_parallel
#  SOURCES  test_solid_mechanics_model_pbc_parallel.cc
#  DEPENDS test_pbc_parallel_mesh
#  FILES_TO_COPY  material.dat
#  DIRECTORIES_TO_CREATE paraview
#  PACKAGE parallel
#  )

#===============================================================================
add_mesh(test_cube3d_mesh1 cube.geo 3 1 OUTPUT cube1.msh)
add_mesh(test_cube3d_mesh2 cube.geo 3 2 OUTPUT cube2.msh)
add_mesh(test_cube3d_mesh_structured cube_structured.geo 3 1 OUTPUT cube_structured.msh)

register_test(test_solid_mechanics_model_cube3d
  SOURCES test_solid_mechanics_model_cube3d.cc
  DEPENDS test_cube3d_mesh1
  FILES_TO_COPY  material.dat
  DIRECTORIES_TO_CREATE paraview
  PACKAGE core
  )

register_test(test_solid_mechanics_model_cube3d_tetra10
  SOURCES test_solid_mechanics_model_cube3d_tetra10.cc
  DEPENDS test_cube3d_mesh2
  FILES_TO_COPY  material.dat
  DIRECTORIES_TO_CREATE paraview
  PACKAGE core
  )

# register_test(test_solid_mechanics_model_cube3d_pbc
#   SOURCES test_solid_mechanics_model_cube3d_pbc.cc
#   DEPENDS test_cube3d_mesh_structured
#   FILES_TO_COPY  material.dat
#   DIRECTORIES_TO_CREATE paraview
#   PACKAGE core
#   )

#add_mesh(test_solid_mechanics_model_boundary_condition_mesh cube_physical_names.geo 3 1)

# register_test(test_solid_mechanics_model_boundary_condition
#   SOURCES test_solid_mechanics_model_boundary_condition.cc
#   DEPENDS test_cube3d_mesh1
#   FILES_TO_COPY material.dat
#   PACKAGE core
#   )

#===============================================================================
add_mesh(test_cube3d_two_mat_mesh cube_two_materials.geo 3 1)

register_test(test_solid_mechanics_model_reassign_material
  SOURCES test_solid_mechanics_model_reassign_material.cc
  DEPENDS test_cube3d_two_mat_mesh
  FILES_TO_COPY two_materials.dat
  PACKAGE implicit
  )

#===============================================================================
register_test(test_solid_mechanics_model_material_eigenstrain
  SOURCES test_solid_mechanics_model_material_eigenstrain.cc
  FILES_TO_COPY cube_3d_tet_4.msh; material_elastic_plane_strain.dat
  PACKAGE core
  )

#===============================================================================
register_test(test_material_selector
  SOURCES test_material_selector.cc
  FILES_TO_COPY material_selector.dat material_selector.msh
  PACKAGE core
  )
