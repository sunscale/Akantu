#===============================================================================
# @file   CMakeLists.txt
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Thu Aug 09 2012
# @date last modification: Wed Feb 03 2016
#
# @brief  configuration for viscoelastic materials tests
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

add_mesh(test_material_viscoelastic_maxwell_mesh
  test_material_viscoelastic_maxwell.geo 2 1)

register_test(test_material_viscoelasti_maxwell_relaxation
  SOURCES test_material_viscoelastic_maxwell_relaxation.cc
  DEPENDS test_material_viscoelastic_maxwell_mesh
  FILES_TO_COPY material_viscoelastic_maxwell.dat
  PACKAGE core
  UNSTABLE
  )
