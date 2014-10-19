#===============================================================================
# @file   CMakeLists.txt
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Nov 28 2011
# @date last modification: Tue Sep 16 2014
#
# @brief  CMake file for the library
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
#===============================================================================

#===============================================================================
# Package Management
#===============================================================================
generate_source_list_from_packages(${CMAKE_CURRENT_BINARY_DIR}
  AKANTU_LIBRARY_SRCS
  AKANTU_LIBRARY_INLINE_HDRS
  AKANTU_LIBRARY_HDRS
  AKANTU_LIBRARY_INCLUDE_DIRS)

#===========================================================================
# header for blas/lapack (any other fortran libraries)
#===========================================================================
if(AKANTU_LAPACK OR AKANTU_BLAS)
  include(FortranCInterface)
  FortranCInterface_HEADER(
    ${CMAKE_CURRENT_BINARY_DIR}/aka_fortran_mangling.hh
    MACRO_NAMESPACE "AKA_FC_")
  mark_as_advanced(CDEFS)
  list(APPEND AKANTU_LIBRARY_HDRS ${CMAKE_CURRENT_BINARY_DIR}/aka_fortran_mangling.hh)
endif()

configure_file(common/aka_config.hh.in "${CMAKE_CURRENT_BINARY_DIR}/aka_config.hh" @ONLY)
list(APPEND AKANTU_LIBRARY_HDRS ${CMAKE_CURRENT_BINARY_DIR}/aka_config.hh)

#===========================================================================
# header precompilation
#===========================================================================
set(AKANTU_COMMON_HDR_TO_PRECOMPILE
  common/aka_vector.hh
  common/aka_math.hh
  common/aka_types.hh
  fem/element_class.hh
  )

set(AKANTU_PRECOMPILE_HDR_INCLUDE_DIRS
  ${CMAKE_CURRENT_BINARY_DIR}/common/
  ${CMAKE_CURRENT_BINARY_DIR}/fem/
  )

set(AKANTU_INCLUDE_DIRS
  ${CMAKE_CURRENT_BINARY_DIR} ${AKANTU_LIBRARY_INCLUDE_DIRS} ${AKANTU_PRECOMPILE_HDR_INCLUDE_DIRS}
  CACHE INTERNAL "Internal include directories to link with Akantu as a subproject")

include_directories(${AKANTU_INCLUDE_DIRS} ${AKANTU_EXTERNAL_LIB_INCLUDE_DIR})

function(generate_material_list)
  message(STATUS "Determining the list of recognized materials...")

  set(_include_dirs ${AKANTU_INCLUDE_DIRS} ${AKANTU_EXTERNAL_LIB_INCLUDE_DIR})

  try_run(_material_list_run _material_list_compile
    ${CMAKE_BINARY_DIR}
    ${PROJECT_SOURCE_DIR}/cmake/material_lister.cc
    CMAKE_FLAGS "-DINCLUDE_DIRECTORIES:STRING=${_include_dirs}"
    COMPILE_DEFINITIONS "-DAKANTU_CMAKE_LIST_MATERIALS"
    COMPILE_OUTPUT_VARIABLE _compile_results
    RUN_OUTPUT_VARIABLE _result_material_list)

  if(_material_list_compile AND "${_material_list_run}" EQUAL 0)
    message(STATUS "Materials included in Akantu:")
    string(REPLACE "\n" ";" _material_list "${_result_material_list}")
    foreach(_mat ${_material_list})
      string(REPLACE ":" ";" _mat_key "${_mat}")
      list(GET _mat_key 0 _key)
      list(GET _mat_key 1 _class)
      list(LENGTH _mat_key _l)

      if("${_l}" GREATER 2)
	list(REMOVE_AT _mat_key 0 1)
	set(_opt " -- options: [")
	foreach(_o ${_mat_key})
	  set(_opt "${_opt} ${_o}")
	endforeach()
	set(_opt "${_opt} ]")
      else()
	set(_opt "")
      endif()

      message(STATUS "   ${_class} -- key: ${_key}${_opt}")
    endforeach()
  else()
    message(STATUS "Could not determine the list of materials.")
    message("${_compile_results}")
  endif()
endfunction()

generate_material_list()

if(CMAKE_COMPILER_IS_GNUCXX)
  include(PCHgcc)
  foreach(_header ${AKANTU_COMMON_HDR_TO_PRECOMPILE})
    add_pch_rule(${_header} AKANTU_LIBRARY_SRCS)
  endforeach()
elseif(CMAKE_COMPILER_IS_GNUCXX)
endif()

#===============================================================================
# Library generation
#===============================================================================
add_library(akantu ${AKANTU_LIBRARY_SRCS})
target_link_libraries(akantu ${AKANTU_EXTERNAL_LIBRARIES})

if(AKANTU_EXTRA_TARGET_DEPENDENCIES)
  add_dependencies(akantu ${AKANTU_EXTRA_TARGET_DEPENDENCIES})
endif()

list(APPEND AKANTU_PUBLIC_HDRS ${AKANTU_LIBRARY_HDRS} ${AKANTU_LIBRARY_INLINE_HDRS})
set_target_properties(akantu PROPERTIES ${AKANTU_LIBRARY_PROPERTIES})
set_target_properties(akantu PROPERTIES PUBLIC_HEADER "${AKANTU_PUBLIC_HDRS}")

list(APPEND AKANTU_EXPORT_LIST akantu)

install(TARGETS akantu
  EXPORT ${AKANTU_TARGETS_EXPORT}
  LIBRARY DESTINATION lib COMPONENT lib
  ARCHIVE DESTINATION lib COMPONENT lib
  PUBLIC_HEADER DESTINATION include/akantu/ COMPONENT dev
  )

if("${AKANTU_TARGETS_EXPORT}" STREQUAL "AkantuLibraryDepends")
  install(EXPORT AkantuLibraryDepends DESTINATION  lib/akantu
    COMPONENT dev)
endif()

#Export for build tree
export(TARGETS ${AKANTU_EXPORT_LIST}
  FILE "${CMAKE_BINARY_DIR}/AkantuLibraryDepends.cmake")
export(PACKAGE Akantu)
