#===============================================================================
# @file   CMakeLists.txt
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
#
# @date creation: Fri Oct 22 2010
# @date last modification: Mon Jan 29 2018
#
# @brief  configuration for materials tests
#
# @section LICENSE
#
# Copyright (©) 2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
# @section DESCRIPTION
#
#===============================================================================

#add_mesh(test_local_material_barre_trou_mesh barre_trou.geo 2 2)
add_mesh(test_local_material_barre_trou_mesh mesh_section_gap.geo 2 2)
register_test(test_local_material
  SOURCES  test_local_material.cc local_material_damage.cc
  EXTRA_FILES local_material_damage.hh local_material_damage_inline_impl.cc
  DEPENDS test_local_material_barre_trou_mesh
  FILES_TO_COPY material.dat
  DIRECTORIES_TO_CREATE paraview
  PACKAGE core
  )

# ==============================================================================
add_mesh(test_interpolate_stress_mesh interpolation.geo 3 2)
register_test(test_interpolate_stress test_interpolate_stress.cc
  FILES_TO_COPY material_interpolate.dat
  DEPENDS test_interpolate_stress_mesh
  DIRECTORIES_TO_CREATE paraview
  PACKAGE lapack core
  )

#===============================================================================
add_mesh(test_material_orthotropic_square_mesh square.geo 2 1)
register_test(test_material_orthotropic
  SOURCES test_material_orthotropic.cc
  DEPENDS test_material_orthotropic_square_mesh
  FILES_TO_COPY orthotropic.dat
  DIRECTORIES_TO_CREATE paraview
  PACKAGE core lapack
  )

#===============================================================================
register_test(test_material_mazars
  SOURCES test_material_mazars.cc
  FILES_TO_COPY material_mazars.dat
  DIRECTORIES_TO_CREATE paraview
  PACKAGE core lapack
  UNSTABLE
  )

# ==============================================================================
add_akantu_test(test_material_viscoelastic "test the visco elastic materials")
add_akantu_test(test_material_non_local "test the non-local materials")
add_akantu_test(test_material_elasto_plastic_linear_isotropic_hardening
  "test the elasto plastic with linear isotropic hardening materials")
add_akantu_test(test_material_viscoelastic_maxwell "test the viscoelastic maxwell material")

# ==============================================================================
add_mesh(test_multi_material_elastic_mesh test_multi_material_elastic.geo 2 1)
register_test(test_multi_material_elastic
  SOURCES test_multi_material_elastic.cc
  FILES_TO_COPY test_multi_material_elastic.dat
  DEPENDS test_multi_material_elastic_mesh
  PACKAGE implicit)


# ==============================================================================
# Material unit tests
# ==============================================================================
add_mesh(test_finite_deformation_mesh cube.geo 3 1)
register_gtest_sources(test_finite_deformation
  SOURCES test_finite_deformation.cc
  FILES_TO_COPY material_finite_deformation.dat
  DEPENDS test_finite_deformation_mesh
  PACKAGE implicit)

register_gtest_sources(SOURCES test_elastic_materials.cc PACKAGE core)
register_gtest_sources(SOURCES test_finite_def_materials.cc PACKAGE core)
register_gtest_sources(SOURCES test_damage_materials.cc PACKAGE core python_interface
  LINK_LIBRARIES pyakantu
  FILES_TO_COPY py_mazars.py)
register_gtest_sources(SOURCES test_plastic_materials.cc PACKAGE core)
register_gtest_sources(SOURCES test_material_thermal.cc PACKAGE core)

register_gtest_test(test_material)
