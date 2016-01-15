#===============================================================================
# @file   FindMumps.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Dec 13 2010
# @date last modification: Tue Sep 09 2014
#
# @brief  The find_package file for the Mumps solver
#
# @section LICENSE
#
# Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

set(_MUMPS_COMPONENTS "sequential" "parallel")

if(NOT Mumps_FIND_COMPONENTS)
  set(Mumps_FIND_COMPONENTS "parallel")
endif()
#===============================================================================
enable_language(Fortran)

if("${Mumps_FIND_COMPONENTS}" STREQUAL "sequential")
  set(MUMPS_PREFIX _seq)
else()
  unset(MUMPS_PREFIX)
endif()

find_path(MUMPS_INCLUDE_DIR dmumps_c.h
  HINTS ${MUMPS_DIR}
  PATH_SUFFIXES include
  )

find_library(MUMPS_LIBRARY_COMMON NAMES mumps_common${MUMPS_PREFIX}
   HINTS ${MUMPS_DIR}
   PATH_SUFFIXES lib
   )

find_library(MUMPS_LIBRARY_PORD NAMES pord${MUMPS_PREFIX}
   HINTS ${MUMPS_DIR}
   PATH_SUFFIXES lib
   )

foreach(_precision s d c z)
  string(TOUPPER "${_precision}" _u_precision)
  find_library(MUMPS_LIBRARY_${_u_precision}MUMPS NAMES ${_precision}mumps${MUMPS_PREFIX}
    HINTS ${MUMPS_DIR} /usr
    PATH_SUFFIXES lib
    )
  mark_as_advanced(MUMPS_LIBRARY_${_u_precision}MUMPS)
endforeach()

mark_as_advanced(
  MUMPS_LIBRARY_COMMON
  MUMPS_LIBRARY_PORD
  MUMPS_INCLUDE_DIR)

#===============================================================================
if(NOT MUMPS_FOUND)
  set(MUMPS_DIR "" CACHE PATH "Prefix of MUMPS library.")
  mark_as_advanced(MUMPS_DIR)
endif()
#===============================================================================
include(FindPackageHandleStandardArgs)
if(CMAKE_VERSION VERSION_GREATER 2.8.12)
  if(MUMPS_INCLUDE_DIR)
    file(STRINGS ${MUMPS_INCLUDE_DIR}/dmumps_c.h _versions
      REGEX "^#define MUMPS_VERSION .*")
    foreach(_ver ${_versions})
      string(REGEX MATCH "MUMPS_VERSION *\"([0-9.]+)\"" _tmp "${_ver}")
      set(_mumps_VERSION ${CMAKE_MATCH_1})
    endforeach()
    set(MUMPS_VERSION "${_mumps_VERSION}" CACHE INTERNAL "")
  endif()

  find_package_handle_standard_args(Mumps
    REQUIRED_VARS
      MUMPS_LIBRARY_DMUMPS
      MUMPS_LIBRARY_COMMON
      MUMPS_LIBRARY_PORD
      MUMPS_INCLUDE_DIR
    VERSION_VAR
      MUMPS_VERSION)
else()
  find_package_handle_standard_args(Mumps DEFAULT_MSG
    MUMPS_LIBRARIES MUMPS_INCLUDE_DIR)
endif()


if (MUMPS_FOUND AND NOT TARGET MUMPS::common)
  set(MUMPS_LIBRARIES_ALL ${MUMPS_LIBRARY_DMUMPS})

  if(MUMPS_LIBRARY_COMMON MATCHES ".*mumps_common.*${CMAKE_STATIC_LIBRARY_SUFFIX}")
    # Assuming mumps was compiled as a static library
    set(MUMPS_LIBRARY_TYPE STATIC CACHE INTERNAL "" FORCE)

    set(_extra_dep_list pthread)
    find_package(BLAS REQUIRED)
    list(APPEND _extra_dep_list ${BLAS_LIBRARIES})

    if (CMAKE_Fortran_COMPILER MATCHES ".*gfortran")
      set(_compiler_specific gfortran)
    elseif (CMAKE_Fortran_COMPILER MATCHES ".*ifort")
      set(_compiler_specific ifcore)
    endif()
    list(APPEND _extra_dep_list ${_compiler_specific})

    list(APPEND MUMPS_LIBRARIES_ALL
      ${MUMPS_LIBRARY_COMMON}
      ${MUMPS_LIBRARY_PORD}
      pthread
      ${_compiler_specific}
      )

    if("${Mumps_FIND_COMPONENTS}" STREQUAL "parallel")
      find_package(MPI REQUIRED)
      list(APPEND _extra_dep_list ${MPI_Fortran_LIBRARIES})

      find_package(ScaLAPACK REQUIRED)
      list(APPEND _extra_dep_list ScaLAPACK)

      list(APPEND MUMPS_LIBRARIES_ALL
        ${MPI_Fortran_LIBRARIES}
        ${SCALAPACK_LIBRARIES}
        )
    endif()

    list(APPEND MUMPS_LIBRARIES_ALL
      ${BLAS_LIBRARIES})

    if(_extra_dep_list)
      set(_extra_dep ";${_extra_dep_list}")
    else()
      set(_extra_dep)
    endif()
  else()
    set(MUMPS_LIBRARY_TYPE SHARED CACHE INTERNAL "" FORCE)
  endif()

  add_library(MUMPS::common ${MUMPS_LIBRARY_TYPE} IMPORTED GLOBAL)
  add_library(MUMPS::pord   ${MUMPS_LIBRARY_TYPE} IMPORTED GLOBAL)

  #TODO adapt it for windows and dlls (check FindGSL as an example)
  set_target_properties(MUMPS::pord PROPERTIES
    IMPORTED_LOCATION                 "${MUMPS_LIBRARY_PORD}"
    INTERFACE_INCLUDE_DIRECTORIES     "${MUMPS_INCLUDE_DIR}"
    IMPORTED_LINK_INTERFACE_LANGUAGES "C")
  set_target_properties(MUMPS::common PROPERTIES
    IMPORTED_LOCATION                 "${MUMPS_LIBRARY_COMMON}"
    INTERFACE_INCLUDE_DIRECTORIES     "${MUMPS_INCLUDE_DIR}"
    IMPORTED_LINK_INTERFACE_LANGUAGES "C;Fortran"
    INTERFACE_LINK_LIBRARIES          "MUMPS::pord${_extra_dep}")

  foreach(_precision s d c z)
    string(TOUPPER "${_precision}" _u_precision)
    set(_target MUMPS::${_precision}mumps)
    add_library(${_target} ${MUMPS_LIBRARY_TYPE} IMPORTED GLOBAL)
    set_target_properties(${_target} PROPERTIES
      IMPORTED_LOCATION                 "${MUMPS_LIBRARY_${_u_precision}MUMPS}"
      INTERFACE_INCLUDE_DIRECTORIES     "${MUMPS_INCLUDE_DIR}"
      IMPORTED_LINK_INTERFACE_LANGUAGES "C;Fortran"
      INTERFACE_LINK_LIBRARIES          "MUMPS::common")
  endforeach()

  set(MUMPS_LIBRARIES ${MUMPS_LIBRARIES_ALL} CACHE INTERNAL "Libraries for MUMPS" FORCE)
endif()
