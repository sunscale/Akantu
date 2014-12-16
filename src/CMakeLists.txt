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
package_get_source_files(
  AKANTU_LIBRARY_SRCS
  AKANTU_LIBRARY_PUBLIC_HDRS
  AKANTU_LIBRARY_PRIVATE_HDRS
  )

package_get_include_directories(
  AKANTU_LIBRARY_INCLUDE_DIRS
  )

package_get_external_informations(
  AKANTU_EXTERNAL_INCLUDE_DIR
  AKANTU_EXTERNAL_LIBRARIES
)

#===========================================================================
# header for blas/lapack (any other fortran libraries)
#===========================================================================
package_get_name(BLAS _blas_pkg_name)
package_get_name(LAPACK _lapack_pkg_name)

package_is_activated(${_blas_pkg_name} _blas_activated)
package_is_activated(${_lapack_pkg_name} _lapack_activated)

if(_blas_activated OR _lapack_activated)
  if(CMAKE_Fortran_COMPILER)
    # ugly hack
    set(CMAKE_Fortran_COMPILER_LOADED TRUE)
  endif()

  include(FortranCInterface)
  FortranCInterface_HEADER(
    ${CMAKE_CURRENT_BINARY_DIR}/aka_fortran_mangling.hh
    MACRO_NAMESPACE "AKA_FC_")
  mark_as_advanced(CDEFS)
  list(APPEND AKANTU_LIBRARY_PUBLIC_HDRS ${CMAKE_CURRENT_BINARY_DIR}/aka_fortran_mangling.hh)
endif()

configure_file(common/aka_config.hh.in "${CMAKE_CURRENT_BINARY_DIR}/aka_config.hh" @ONLY)
list(APPEND AKANTU_LIBRARY_PUBLIC_HDRS ${CMAKE_CURRENT_BINARY_DIR}/aka_config.hh)
list(APPEND AKANTU_LIBRARY_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR})
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

package_get_all_extra_dependency(_extra_target_dependencies)
if(_extra_target_dependencies)
  # This only adding todo: find a solution for when a dependency was add the is removed...
  add_dependencies(akantu ${_extra_target_dependencies})
endif()

set_target_properties(akantu PROPERTIES ${AKANTU_LIBRARY_PROPERTIES})
set_target_properties(akantu PROPERTIES PUBLIC_HEADER "${AKANTU_LIBRARY_PUBLIC_HDRS}")

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
