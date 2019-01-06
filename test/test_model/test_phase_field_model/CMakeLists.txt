#===============================================================================
# @file   CMakeLists.txt
#
# @author Mohit Pundir <mohit.pundir@epfl.ch>
#
# @date creation: Thu Dec 20 2018
# @date last modification: Thu Dec 20 2018
#
# @brief  test for the phase field model
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
#===============================================================================

add_mesh(test_one_element test_one_element.geo 2 1)
register_test(test_phase_field_model_2d
  SOURCES test_phase_field_model_2d.cc
  DEPENDS test_one_element
  FILES_TO_COPY material.dat
  PACKAGE phase_field
  )

register_test(test_solid_one_element
  SOURCES test_solid_one_element.cc
  DEPENDS test_one_element
  FILES_TO_COPY material_solid.dat
  PACKAGE core
  )

register_test(test_phase_solid_coupling
  SOURCES test_phase_solid_coupling.cc
  DEPENDS test_one_element
  FILES_TO_COPY material_coupling.dat
  PACKAGE phase_field core
  )
