#===============================================================================
# @file   CMakeLists.txt
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @author Clement Roux <clement.roux@epfl.ch>
#
# @date creation: Fri Sep 03 2010
# @date last modification: Wed Feb 03 2016
#
# @brief  configuration for materials tests
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

add_mesh(test_material_non_local_mesh mesh.geo 2 1 OUTPUT mesh.msh)
add_mesh(test_material_damage_non_local_mesh mesh_section_gap.geo 2 1 OUTPUT mesh_section_gap.msh)

register_test(test_material_damage_non_local
  SOURCES test_material_damage_non_local.cc
  DEPENDS test_material_damage_non_local_mesh
  FILES_TO_COPY material_damage_non_local.dat
  DIRECTORIES_TO_CREATE paraview
  PACKAGE damage_non_local
  )


register_test(test_material_non_local
  SOURCES test_material_non_local.cc custom_non_local_test_material.cc custom_non_local_test_material.hh
  DEPENDS test_material_non_local_mesh
  FILES_TO_COPY material.dat
  PACKAGE damage_non_local
  )
