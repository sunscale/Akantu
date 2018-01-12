#===============================================================================
# @file   CMakeLists.txt
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
#
# @date creation: Fri Sep 03 2010
# @date last modification: Mon Dec 07 2015
#
# @brief  configuration for FEM tests
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
package_get_element_types(core _types)

foreach(_type ${_types})
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${_type}.msh)
    list(APPEND _meshes ${_type}.msh)
    #register_fem_test(fe_engine_precomputation ${_type})
  else()
    if(NOT ${_type} STREQUAL _point_1)
      message("The mesh ${_type}.msh is missing, the fe_engine test cannot be activated without it")
    endif()
  endif()
endforeach()


register_test(test_fe_engine_precomupation_bernoulli_beam_2
  PACKAGE structural_mechanics
  SOURCES test_fe_engine_precomputation_bernoulli_2.cc)
#add_mesh(test_fem_circle_1_mesh circle.geo 2 1 OUTPUT circle1.msh)
#add_mesh(test_fem_circle_2_mesh circle.geo 2 2 OUTPUT circle2.msh)

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

package_is_activated(pybind11 _pybind11_act)
if(_pybind11_act)
  pybind11_add_module(aka_test MODULE pybind11_akantu.cc)
endif()

register_gtest_sources(
  SOURCES test_fe_engine_precomputation.cc
  PACKAGE core pybind11
  DEPENDS aka_test
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
