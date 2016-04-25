#===============================================================================
# @file   FindScotch.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Oct 24 2014
# @date last modification: Wed Jan 13 2016
#
# @brief  The find_package file for Scotch
#
# @section LICENSE
#
# Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
# (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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

set(_SCOTCH_COMPONENTS "metis" "parmetis" "esmumps" "ptscotch")

if(NOT Scotch_FIND_COMPONENTS)
  set(Scotch_FIND_COMPONENTS)
endif()

find_path(SCOTCH_INCLUDE_DIR scotch.h  PATHS "${SCOTCH_DIR}" ENV SCOTCH_DIR
  PATH_SUFFIXES include include/scotch
  )


find_library(SCOTCH_LIBRARY scotch PATHS "${SCOTCH_DIR}" ENV SCOTCH_DIR PATH_SUFFIXES lib)

set(_scotch_test_dir "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}")
file(WRITE "${_scotch_test_dir}/scotch_test_code.c"
  "#include <stdio.h>
#include <stdint.h>
#include <scotch.h>

int main() {
  SCOTCH_Graph graph;
  SCOTCH_graphInit(&graph);
  return 0;
}
")

#===============================================================================
include(FindPackageHandleStandardArgs)
if(CMAKE_VERSION VERSION_GREATER 2.8.12)
  if(SCOTCH_INCLUDE_DIR)
    file(STRINGS ${SCOTCH_INCLUDE_DIR}/scotch.h _versions
      REGEX "^#define\ +SCOTCH_(VERSION|RELEASE|PATCHLEVEL) .*")
    foreach(_ver ${_versions})
      string(REGEX MATCH "SCOTCH_(VERSION|RELEASE|PATCHLEVEL) *([0-9.]+)" _tmp "${_ver}")
      set(_scotch_${CMAKE_MATCH_1} ${CMAKE_MATCH_2})
    endforeach()
    set(SCOTCH_VERSION "${_scotch_VERSION}.${_scotch_RELEASE}.${_scotch_PATCHLEVEL}" CACHE INTERNAL "")
  endif()
  find_package_handle_standard_args(Scotch
    REQUIRED_VARS SCOTCH_LIBRARY SCOTCH_INCLUDE_DIR
    VERSION_VAR SCOTCH_VERSION)
else()
  find_package_handle_standard_args(Scotch DEFAULT_MSG
    SCOTCH_LIBRARY SCOTCH_INCLUDE_DIR)
endif()

set(SCOTCH_LIBRARIES_ALL ${SCOTCH_LIBRARY})

try_compile(_scotch_compiles "${_scotch_test_dir}" SOURCES "${_scotch_test_dir}/scotch_test_code.c"
  CMAKE_FLAGS "-DINCLUDE_DIRECTORIES:STRING=${SCOTCH_INCLUDE_DIR}"
  LINK_LIBRARIES ${SCOTCH_LIBRARY}
  OUTPUT_VARIABLE _out)

get_filename_component(_scotch_hint "${SCOTCH_LIBRARY}" DIRECTORY)

if(SCOTCH_LIBRARY MATCHES ".*scotch.*${CMAKE_STATIC_LIBRARY_SUFFIX}")
  # Assuming scotch was compiled as a static library
  set(SCOTCH_LIBRARY_TYPE STATIC CACHE INTERNAL "" FORCE)
else()
  set(SCOTCH_LIBRARY_TYPE SHARED CACHE INTERNAL "" FORCE)
endif()

if(NOT _scotch_compiles)
  if(_out MATCHES "SCOTCH_errorPrint")
    find_library(SCOTCH_LIBRARY_ERR scotcherr
      HINTS ${_scotch_hint})
    find_library(SCOTCH_LIBRARY_ERREXIT scotcherrexit
      HINTS ${_scotch_hint})

    if(NOT TARGET Scotch::err)
      add_library(Scotch::err ${SCOTCH_LIBRARY_TYPE} IMPORTED GLOBAL)
    endif()
    if(NOT TARGET Scotch::errexit)
      add_library(Scotch::errexit ${SCOTCH_LIBRARY_TYPE} IMPORTED GLOBAL)
    endif()

    set_target_properties(Scotch::errexit PROPERTIES
      IMPORTED_LOCATION                 "${SCOTCH_LIBRARY_ERREXIT}"
      INTERFACE_INCLUDE_DIRECTORIES     "${SCOTCH_INCLUDE_DIR}"
      IMPORTED_LINK_INTERFACE_LANGUAGES "C")

    set_target_properties(Scotch::err PROPERTIES
      IMPORTED_LOCATION                 "${SCOTCH_LIBRARY_ERR}"
      INTERFACE_INCLUDE_DIRECTORIES     "${SCOTCH_INCLUDE_DIR}"
      IMPORTED_LINK_INTERFACE_LANGUAGES "C"
      INTERFACE_LINK_LIBRARIES          "Scotch::errexit")

    mark_as_advanced(SCOTCH_LIBRARY_ERR
      SCOTCH_LIBRARY_ERREXIT)

    list(APPEND SCOTCH_LIBRARIES_ALL ${SCOTCH_LIBRARY_ERR} ${SCOTCH_LIBRARY_ERREXIT})

    set(_scotch_link_lib INTERFACE_LINK_LIBRARIES "Scotch::err")
  endif()
endif()

if(NOT TARGET Scotch::scotch)
  add_library(Scotch::scotch ${SCOTCH_LIBRARY_TYPE} IMPORTED GLOBAL)
endif()
set_target_properties(Scotch::scotch PROPERTIES
  IMPORTED_LOCATION                 "${SCOTCH_LIBRARY}"
  INTERFACE_INCLUDE_DIRECTORIES     "${SCOTCH_INCLUDE_DIR}"
  IMPORTED_LINK_INTERFACE_LANGUAGES "C"
  ${_scotch_link_lib})

set(SCOTCH_LIBRARIES ${SCOTCH_LIBRARIES_ALL} CACHE INTERNAL "Libraries for Scotch" FORCE)

mark_as_advanced(SCOTCH_LIBRARY
  SCOTCH_INCLUDE_DIR
  SCOTCH_LIBRARIES)


if("${Scotch_FIND_COMPONENTS}" MATCHES "esmumps")
  find_library(SCOTCH_LIBRARY_ESMUMPS esmumps  HINTS ${_scotch_hint})

  if(NOT TARGET Scotch::esmumps)
    add_library(Scotch::esmumps ${SCOTCH_LIBRARY_TYPE} IMPORTED GLOBAL)
  endif()
  set_target_properties(Scotch::esmumps PROPERTIES
    IMPORTED_LOCATION                 "${SCOTCH_LIBRARY_ESMUMPS}"
    INTERFACE_INCLUDE_DIRECTORIES     "${SCOTCH_INCLUDE_DIR}"
    IMPORTED_LINK_INTERFACE_LANGUAGES "C")


  mark_as_advanced(SCOTCH_LIBRARY_ESMUMPS)
endif()

if("${Scotch_FIND_COMPONENTS}" MATCHES "metis")
  find_library(SCOTCH_LIBRARY_METIS scotchmetis HINTS ${_scotch_hint})

  if(NOT TARGET Scotch::metis)
    add_library(Scotch::metis ${SCOTCH_LIBRARY_TYPE} IMPORTED GLOBAL)
  endif()
  set_target_properties(Scotch::metis PROPERTIES
    IMPORTED_LOCATION                 "${SCOTCH_LIBRARY_METIS}"
    INTERFACE_INCLUDE_DIRECTORIES     "${SCOTCH_INCLUDE_DIR}"
    IMPORTED_LINK_INTERFACE_LANGUAGES "C")

  mark_as_advanced(SCOTCH_LIBRARY_METIS)
endif()

if("${Scotch_FIND_COMPONENTS}" MATCHES "parmetis")
  find_library(SCOTCH_LIBRARY_PARMETIS scotchparmetis HINTS ${_scotch_hint})

  if(NOT TARGET Scotch::parmetis)
    add_library(Scotch::parmetis ${SCOTCH_LIBRARY_TYPE} IMPORTED GLOBAL)
  endif()
  set_target_properties(Scotch::parmetis PROPERTIES
    IMPORTED_LOCATION                 "${SCOTCH_LIBRARY_PARMETIS}"
    INTERFACE_INCLUDE_DIRECTORIES     "${SCOTCH_INCLUDE_DIR}"
    IMPORTED_LINK_INTERFACE_LANGUAGES "C"
    INTERFACE_INCLUDE_DIRECTORIES     "Scotch::metis")
  mark_as_advanced(SCOTCH_LIBRARY_PARMETIS)
endif()

#
##===============================================================================
if("${Scotch_FIND_COMPONENTS}" MATCHES "ptscotch")
  file(WRITE "${_scotch_test_dir}/ptscotch_test_code.c"
    "#include <stdio.h>
#include <stdint.h>
#include <mpi.h>
#include <ptscotch.h>

int main() {
  SCOTCH_Dgraph graph;
  SCOTCH_dgraphInit(&graph, MPI_COMM_WORLD);
  return 0;
}
")

  find_package(MPI REQUIRED)

  find_library(SCOTCH_LIBRARY_PTSCOTCH ptscotch HINTS ${_scotch_hint})

  try_compile(_scotch_compiles "${_scotch_test_dir}" SOURCES "${_scotch_test_dir}/ptscotch_test_code.c"
    CMAKE_FLAGS "-DINCLUDE_DIRECTORIES:STRING=${SCOTCH_INCLUDE_DIR};${MPI_C_INCLUDE_PATH}"
    LINK_LIBRARIES ${SCOTCH_LIBRARY_PTSCOTCH} ${MPI_C_LIBRARIES}
    OUTPUT_VARIABLE _out)

  if(NOT _scotch_compiles)
    if(_out MATCHES "SCOTCH_archExit")
      set(_scotch_link_lib INTERFACE_LINK_LIBRARIES "Scotch::scotch")
    endif()
  endif()

  if(NOT TARGET Scotch::ptscotch)
    add_library(Scotch::ptscotch ${SCOTCH_LIBRARY_TYPE} IMPORTED GLOBAL)
  endif()
  set_target_properties(Scotch::ptscotch PROPERTIES
    IMPORTED_LOCATION                 "${SCOTCH_LIBRARY_PTSCOTCH}"
    INTERFACE_INCLUDE_DIRECTORIES     "${SCOTCH_INCLUDE_DIR}"
    IMPORTED_LINK_INTERFACE_LANGUAGES "C"
    ${_scotch_link_lib})

  set(PTSCOTCH_LIBRARIES ${SCOTCH_LIBRARY_PTSCOTCH} ${SCOTCH_LIBRARIES} CACHE INTERNAL "Libraries for PT-Scotch" FORCE)

  mark_as_advanced(SCOTCH_LIBRARY_PTSCOTCH
    PTSCOTCH_LIBRARIES)

  if("${Scotch_FIND_COMPONENTS}" MATCHES "esmumps")
    find_library(SCOTCH_LIBRARY_PTESMUMPS ptesmumps
      HINTS ${_scotch_hint} PATH_SUFFIXES lib .)

    if(NOT TARGET Scotch::ptesmumps)
      add_library(Scotch::ptesmumps ${SCOTCH_LIBRARY_TYPE} IMPORTED GLOBAL)
    endif()
    set_target_properties(Scotch::ptesmumps PROPERTIES
      IMPORTED_LOCATION                 "${SCOTCH_LIBRARY_ESMUMPS}"
      INTERFACE_INCLUDE_DIRECTORIES     "${SCOTCH_INCLUDE_DIR}"
      IMPORTED_LINK_INTERFACE_LANGUAGES "C")

    mark_as_advanced(SCOTCH_LIBRARY_PTESMUMPS)
  endif()
endif()
