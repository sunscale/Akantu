#===============================================================================
# @file   CMakeLists.txt
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Jun 14 2010
# @date last modification: Wed Jan 20 2016
#
# @brief  CMake file for the library
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
#===============================================================================

#===============================================================================
# Package Management
#===============================================================================
package_get_all_source_files(
  AKANTU_LIBRARY_SRCS
  AKANTU_LIBRARY_PUBLIC_HDRS
  AKANTU_LIBRARY_PRIVATE_HDRS
  )

package_get_all_include_directories(
  AKANTU_LIBRARY_INCLUDE_DIRS
  )

package_get_all_external_informations(
  AKANTU_EXTERNAL_INCLUDE_DIR
  AKANTU_EXTERNAL_LIBRARIES
)

package_get_all_compilation_flags(CXX _cxx_flags)
set(AKANTU_EXTRA_CXX_FLAGS "${_cxx_flags}" CACHE STRING "Extra flags defined by loaded packages" FORCE)
mark_as_advanced(AKANTU_EXTRA_CXX_FLAGS)

#===========================================================================
# header for blas/lapack (any other fortran libraries)
#===========================================================================
package_is_activated(BLAS _blas_activated)
package_is_activated(LAPACK _lapack_activated)

if(_blas_activated OR _lapack_activated)
  if(CMAKE_Fortran_COMPILER)
    # ugly hack
    set(CMAKE_Fortran_COMPILER_LOADED TRUE)
  endif()

  include(FortranCInterface)
  FortranCInterface_HEADER(
    "${CMAKE_CURRENT_BINARY_DIR}/aka_fortran_mangling.hh"
    MACRO_NAMESPACE "AKA_FC_")
  mark_as_advanced(CDEFS)
  list(APPEND AKANTU_LIBRARY_PUBLIC_HDRS
    "${CMAKE_CURRENT_BINARY_DIR}/aka_fortran_mangling.hh")
endif()


#===========================================================================
# configurations
#===========================================================================
package_get_all_material_includes(AKANTU_MATERIAL_INCLUDES)
package_get_all_material_lists(AKANTU_MATERIAL_LISTS)
configure_file(model/solid_mechanics/material_list.hh.in
  "${CMAKE_CURRENT_BINARY_DIR}/material_list.hh" @ONLY)

package_get_element_lists()
configure_file(common/aka_element_classes_info.hh.in
  "${CMAKE_CURRENT_BINARY_DIR}/aka_element_classes_info.hh" @ONLY)

configure_file(common/aka_config.hh.in
  "${CMAKE_CURRENT_BINARY_DIR}/aka_config.hh" @ONLY)

list(APPEND AKANTU_LIBRARY_PUBLIC_HDRS
  "${CMAKE_CURRENT_BINARY_DIR}/aka_config.hh"
  "${CMAKE_CURRENT_BINARY_DIR}/aka_element_classes_info.hh"
  "${CMAKE_CURRENT_BINARY_DIR}/material_list.hh")
list(APPEND AKANTU_LIBRARY_INCLUDE_DIRS "${CMAKE_CURRENT_BINARY_DIR}")

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

include_directories(${AKANTU_INCLUDE_DIRS} ${AKANTU_EXTERNAL_INCLUDE_DIR})

if(CMAKE_COMPILER_IS_GNUCXX)
  include(PCHgcc)
  foreach(_header ${AKANTU_COMMON_HDR_TO_PRECOMPILE})
    add_pch_rule(${_header} AKANTU_LIBRARY_SRCS)
  endforeach()
elseif(CMAKE_COMPILER_IS_GNUCXX)
endif()

#===============================================================================
# Debug infos
#===============================================================================
set(AKANTU_GDB_DIR ${PROJECT_SOURCE_DIR}/cmake)
if(UNIX)
  string(TOUPPER "${CMAKE_BUILD_TYPE}" _u_build_type)
  if(_u_build_type STREQUAL "DEBUG" OR _u_build_type STREQUAL "RELWITHDEBINFO")
    configure_file(${PROJECT_SOURCE_DIR}/cmake/libakantu-gdb.py.in
      "${PROJECT_BINARY_DIR}/libakantu-gdb.py"
      @ONLY)
    configure_file(${PROJECT_SOURCE_DIR}/cmake/akantu-debug.cc.in
      "${PROJECT_BINARY_DIR}/akantu-debug.cc" @ONLY)

    list(APPEND AKANTU_LIBRARY_SRCS ${PROJECT_BINARY_DIR}/akantu-debug.cc)
  endif()
else()
  find_program(GDB_EXECUTABLE gdb)

  if(GDB_EXECUTABLE)
    execute_process(COMMAND
      ${GDB_EXECUTABLE} --batch -x "${PROJECT_SOURCE_DIR}/cmake/gdb_python_path"
      OUTPUT_VARIABLE AKANTU_PYTHON_GDB_DIR
      ERROR_QUIET
      RESULT_VARIABLE _res)

    if(_res EQUAL 0 AND UNIX)
      set(GDB_USER_CONFIG $ENV{HOME}/.gdb/auto-load)
      file(MAKE_DIRECTORY ${GDB_USER_CONFIG})

      configure_file(${PROJECT_SOURCE_DIR}/cmake/libakantu-gdb.py.in
        "${GDB_USER_CONFIG}/${CMAKE_SHARED_LIBRARY_PREFIX}akantu${CMAKE_SHARED_LIBRARY_SUFFIX}.${AKANTU_VERSION}-gdb.py"
        @ONLY)
    endif()
  endif()
endif()

#===============================================================================
# Library generation
#===============================================================================
add_library(akantu ${AKANTU_LIBRARY_SRCS})

# the repetition is ugly but works for static libraries
target_link_libraries(akantu  ${AKANTU_EXTERNAL_LIBRARIES} ${AKANTU_EXTERNAL_LIBRARIES})

set_target_properties(akantu
  PROPERTIES
    ${AKANTU_LIBRARY_PROPERTIES} # this contains the version
    PUBLIC_HEADER "${AKANTU_LIBRARY_PUBLIC_HDRS}"
    COMPILE_FLAGS "${_cxx_flags}"
    LINK_INTERFACE_MULTIPLICITY 2)

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.1)
  package_get_all_features_public(_public_features)
  if(_public_features)
    target_compile_features(akantu  PUBLIC ${_public_features})
  endif()

  package_get_all_features_public(_private_features)
  if(_private_features)
    target_compile_features(akantu PRIVATE ${_private_features})
  endif()
endif()

package_get_all_extra_dependencies(_extra_target_dependencies)
if(_extra_target_dependencies)
  # This only adding todo: find a solution for when a dependency was add the is removed...
  add_dependencies(akantu ${_extra_target_dependencies})
endif()


list(APPEND AKANTU_EXPORT_LIST akantu)

# TODO separate public from private headers
install(TARGETS akantu
  EXPORT ${AKANTU_TARGETS_EXPORT}
  LIBRARY DESTINATION lib COMPONENT lib
  ARCHIVE DESTINATION lib COMPONENT lib
  RUNTIME DESTINATION bin COMPONENT bin
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

# print out the list of materials
generate_material_list()
