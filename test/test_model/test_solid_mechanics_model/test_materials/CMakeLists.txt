#===============================================================================
# @file   CMakeLists.txt
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
#
# @date creation: Fri Nov 26 2010
# @date last modification: Fri Feb 14 2014
#
# @brief  configuration for materials tests
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

add_mesh(test_local_material_barre_trou_mesh barre_trou.geo 2 2)
register_test(test_local_material
  SOURCES  test_local_material.cc local_material_damage.cc
  EXTRA_FILES local_material_damage.hh local_material_damage_inline_impl.cc
  DEPENDENCIES test_local_material_barre_trou_mesh
  FILES_TO_COPY material.dat
  DIRECTORIES_TO_CREATE paraview
  )

add_mesh(test_material_thermal_mesh square.geo 2 1)
register_test(test_material_thermal
  SOURCES test_material_thermal.cc
  DEPENDENCIES test_material_thermal_mesh
  FILES_TO_COPY material_thermal.dat
)
#===============================================================================
add_mesh(test_interpolate_stress_mesh interpolation.geo 3 2)
register_test(test_interpolate_stress test_interpolate_stress.cc
  FILES_TO_COPY material.dat
  DEPENDENCIES test_interpolate_stress_mesh
  DIRECTORIES_TO_CREATE paraview)
#===============================================================================

add_akantu_test(test_material_plasticity "test material plasticity"
  PACKAGE extra_materials)

add_akantu_test(test_material_viscoelastic "test the visco elastic materials"
  PACKAGE extra_materials)

add_akantu_test(test_material_non_local "test the non-local materials"
  PACKAGE damage_non_local)
