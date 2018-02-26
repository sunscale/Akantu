#===============================================================================
# @file   CMakeLists.txt
#
# @author Tobias Brink <tobias.brink@epfl.ch>
#
# @date creation: Fri Jan 26 2018
# @date last modification: Tue Jan 30 2018
#
# @brief  
#
# @section LICENSE
#
# Copyright (©) 2016-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free  software: you can redistribute it and/or  modify it under the terms  of the  GNU Lesser  General Public  License as published by  the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more details.
#
# You should  have received  a copy  of the GNU  Lesser General  Public License along with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================

#===============================================================================
# Energy tests
#===============================================================================
register_gtest_sources(
  SOURCES test_solid_mechanics_model_linear_elastic_potential_energy.cc
  FILES_TO_COPY test_solid_mechanics_model_linear_elastic_potential_energy_material.dat
  PACKAGE implicit
  )

register_gtest_sources(
  SOURCES test_solid_mechanics_model_kinetic_energy.cc
  FILES_TO_COPY test_solid_mechanics_model_kinetic_energy_material.dat
  PACKAGE implicit
  )

register_gtest_sources(
  SOURCES test_solid_mechanics_model_work_quasistatic.cc
  FILES_TO_COPY test_solid_mechanics_model_work_material.dat
  PACKAGE implicit
  )

register_gtest_sources(
  SOURCES test_solid_mechanics_model_work_dynamics.cc
  FILES_TO_COPY test_solid_mechanics_model_work_material.dat
  PACKAGE core
  )

register_gtest_test(test_solid_mechanics_model_energies
  FILES_TO_COPY ${PATCH_TEST_MESHES}
  DEPENDS ${PATCH_TEST_BAR_MESHES}
  PARALLEL
  PARALLEL_LEVEL 1 2 # some meshes are too small for more
  )
