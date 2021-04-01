#===============================================================================
# @file   CMakeLists.txt
#
# @author Fabian Barras <fabian.barras@epfl.ch>
#
# @date creation: Fri Sep 03 2010
# @date last modification: Fri Feb 09 2018
#
# @brief  Structural Mechanics Model Tests
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
#===============================================================================

# Adding sources
register_gtest_sources(
  SOURCES test_structural_mechanics_model_bernoulli_beam_2.cc
  PACKAGE implicit structural_mechanics
  )

register_gtest_sources(
  SOURCES test_structural_mechanics_model_bernoulli_beam_3.cc
  PACKAGE implicit structural_mechanics
  )

register_gtest_sources(
  SOURCES test_structural_mechanics_model_discrete_kirchhoff_triangle_18.cc
  PACKAGE implicit structural_mechanics
  )

register_gtest_sources(
  SOURCES test_structural_mechanics_model_bernoulli_beam_dynamics.cc
  PACKAGE implicit structural_mechanics
  )

#===============================================================================

# Adding meshes for element types
package_get_element_types(structural_mechanics _types)

foreach(_type ${_types})
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${_type}.msh)
    list(APPEND _meshes ${_type}.msh)
    #register_fem_test(fe_engine_precomputation ${_type })
  else()
    if(NOT ${_type} STREQUAL _point_1)
      message("The mesh ${_type}.msh is missing, the fe_engine test cannot be activated without it")
    endif()
  endif()
endforeach()

#===============================================================================
# Registering google test
register_gtest_test(test_structural_mechanics
  FILES_TO_COPY ${_meshes}
  )
