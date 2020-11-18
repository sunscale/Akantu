#===============================================================================
# @file   CMakeLists.txt
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Lucas Frerot <lucas.frerot@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Sep 03 2010
# @date last modification: Fri Jan 26 2018
#
# @brief  configuration for FEM tests
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
# You should have received a copy of the GNU Lesser General Public License
# along with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
# @section DESCRIPTION
#
#===============================================================================

#===============================================================================
function(register_fem_test operation type)
  set(_target test_${operation}${type})

  register_test(${_target}
    SOURCES test_${operation}.cc
    FILES_TO_COPY ${type}.msh
    COMPILE_OPTIONS TYPE=${type}
    PACKAGE core
    )
endfunction()

#===============================================================================
macro(register_mesh_types package)
  package_get_element_types(${package} _types)

  foreach(_type ${_types})
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${_type}.msh)
      list(APPEND _meshes ${_type}.msh)
    else()
      if(NOT ${_type} STREQUAL _point_1)
	message("The mesh ${_type}.msh is missing, the fe_engine test cannot be activated without it")
      endif()
    endif()
  endforeach()
endmacro(register_mesh_types)

set(_meshes)
register_mesh_types(core)

package_is_activated(structural_mechanics has_structural_mechanics)
if(has_structural_mechanics)
  register_mesh_types(structural_mechanics)
endif()

# Tests for class MeshData
macro(register_typed_test test_name type value1 value2)
  set(target test_${test_name}_${type})
  register_test(${target}
    SOURCES test_${test_name}.cc
    COMPILE_OPTIONS "TYPE=${type};VALUE1=${value1};VALUE2=${value2}"
    PACKAGE core
    )
endmacro()

register_typed_test(mesh_data string \"5\" \"10\")
register_typed_test(mesh_data UInt 5 10)

add_mesh(test_boundary_msh cube.geo 3 1)
add_mesh(test_boundary_msh_physical_names cube_physical_names.geo 3 1)

register_test(test_mesh_boundary
  SOURCES test_mesh_boundary.cc
  DEPENDS test_boundary_msh test_boundary_msh_physical_names
  PACKAGE core)

register_test(test_facet_element_mapping
  SOURCES test_facet_element_mapping.cc
  DEPENDS test_boundary_msh_physical_names
  PACKAGE core)

register_gtest_sources(
  SOURCES test_fe_engine_precomputation.cc
  PACKAGE core python_interface
  INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/python
  )

register_gtest_sources(
  SOURCES test_fe_engine_precomputation_structural.cc
  PACKAGE structural_mechanics
)


register_gtest_sources(
  SOURCES test_fe_engine_gauss_integration.cc
  PACKAGE core
  )

register_gtest_sources(
  SOURCES test_gradient.cc
  PACKAGE core
  )

register_gtest_sources(
  SOURCES test_integrate.cc
  PACKAGE core
  )

register_gtest_sources(
  SOURCES test_inverse_map.cc
  PACKAGE core
  )

register_gtest_test(test_fe_engine
  FILES_TO_COPY ${_meshes})
