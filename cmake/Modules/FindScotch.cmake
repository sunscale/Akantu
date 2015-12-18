#===============================================================================
# @file   FindScotch.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Wed Sep 01 2010
# @date last modification: Tue Sep 09 2014
#
# @brief  The find_package file for Scotch
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

#===============================================================================
find_library(SCOTCH_LIBRARY scotch
  HINTS ${SCOTCH_DIR}
  PATH_SUFFIXES src/libscotch lib
  )

find_library(SCOTCH_LIBRARY_ERR scotcherr
  HINTS ${SCOTCH_DIR}
  PATH_SUFFIXES src/libscotch lib
  )

find_library(SCOTCH_LIBRARY_ERREXIT scotcherrexit
  HINTS ${SCOTCH_DIR}
  PATH_SUFFIXES src/libscotch lib
  )

find_library(SCOTCH_LIBRARY_ESMUMPS esmumps
  HINTS ${SCOTCH_DIR}
  PATH_SUFFIXES src/libscotch lib
  )

find_library(SCOTCH_LIBRARY_METIS scotchmetis
  PATHS ${SCOTCH_DIR}
  PATH_SUFFIXES src/libscotch lib
  )

find_library(SCOTCH_LIBRARY_PARMETIS scotchparmetis
  PATHS ${SCOTCH_DIR}
  PATH_SUFFIXES src/libscotch lib
  )

find_path(SCOTCH_INCLUDE_DIR scotch.h
  HINTS ${SCOTCH_DIR}
  PATH_SUFFIXES include scotch src/libscotch include/scotch
  )

#===============================================================================
mark_as_advanced(SCOTCH_LIBRARY
  SCOTCH_LIBRARY_ERR
  SCOTCH_LIBRARY_ERREXIT
  SCOTCH_LIBRARY_ESMUMPS
  SCOTCH_LIBRARY_PARMETIS
  SCOTCH_INCLUDE_DIR)

set(SCOTCH_LIBRARIES_ALL ${SCOTCH_LIBRARY} ${SCOTCH_LIBRARY_ERR})

if(SCOTCH_LIBRARY_ESMUMPS)
  set(SCOTCH_LIBRARIES_ALL ${SCOTCH_LIBRARY_ESMUMPS} ${SCOTCH_LIBRARIES_ALL})
endif()

if(SCOTCH_LIBRARY_METIS)
  set(SCOTCH_LIBRARIES_ALL ${SCOTCH_LIBRARY_METIS} ${SCOTCH_LIBRARIES_ALL})
endif()

if(SCOTCH_LIBRARY_PARMETIS)
  set(SCOTCH_LIBRARIES_ALL ${SCOTCH_LIBRARY_PARMETIS} ${SCOTCH_LIBRARIES_ALL})
endif()


set(SCOTCH_LIBRARIES ${SCOTCH_LIBRARIES_ALL} CACHE INTERNAL "Libraries for scotch" FORCE)

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
    set(SCOTCH_VERSION "${_scotch_VERSION}.${_scotch_PATCHLEVEL}" CACHE INTERNAL "")
  endif()
  find_package_handle_standard_args(Scotch
    REQUIRED_VARS SCOTCH_LIBRARIES SCOTCH_INCLUDE_DIR
    VERSION_VAR SCOTCH_VERSION)
else()
  find_package_handle_standard_args(Scotch DEFAULT_MSG
    SCOTCH_LIBRARIES SCOTCH_INCLUDE_DIR)
endif()

#===============================================================================
if(NOT SCOTCH_FOUND)
  set(SCOTCH_DIR "" CACHE PATH "Location of Scotch library.")
  mark_as_advanced(SCOTCH_DIR)
endif()
