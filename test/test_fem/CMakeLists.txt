#===============================================================================
# @file   CMakeLists.txt
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
#
# @date creation: Fri Sep 03 2010
# @date last modification: Fri May 30 2014
#
# @brief  configuration for FEM tests
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

#===============================================================================
macro(register_fem_test operation type)
  set(_target test_${operation}_${_type})
  register_test(${_target}
    SOURCES test_${operation}.cc
    FILES_TO_COPY _${type}.msh
    COMPILE_OPTIONS TYPE=_${_type}
    PACKAGE core
    )
endmacro()
#===============================================================================

set(LIST_TYPES
  segment_2
  segment_3
  triangle_3
  triangle_6
  quadrangle_4
  quadrangle_8
  tetrahedron_4
  tetrahedron_10
  hexahedron_8
  )

foreach(_type ${LIST_TYPES})
  register_fem_test(interpolate ${_type})
  register_fem_test(gradient    ${_type})
  register_fem_test(integrate   ${_type})
  register_fem_test(inverse_map ${_type})
endforeach()


#register_test(test_interpolate_bernoulli_beam_2 test_interpolate_bernoulli_beam_2.cc)
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
  DEPENDENCIES test_boundary_msh test_boundary_msh_physical_names
  PACKAGE core)

register_test(test_facet_element_mapping
  SOURCES test_facet_element_mapping.cc
  DEPENDENCIES test_boundary_msh_physical_names
  PACKAGE core)

register_test(test_igfem_integrate
  SOURCES test_igfem_integrate.cc
  PACKAGE igfem)
