#===============================================================================
# @file   FindPTScotch.cmake
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @date   Tue Aug  25 16:53:57 2010
#
# @brief  The find_package file for PT-Scotch
#
# @section LICENSE
#
# Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#if(PTSCOTCH_DIR)
#  set(PTSCOTCH_LIBRARY "NOTFOUND" CACHE INTERNAL "Cleared" FORCE)
#endif(PTSCOTCH_DIR)

find_library(PTSCOTCH_LIBRARY ptscotch
  PATHS ${PTSCOTCH_DIR}
  PATH_SUFFIXES src/libscotch lib
  )

find_library(PTSCOTCH_LIBRARY_ERR ptscotcherr
  PATHS ${PTSCOTCH_DIR}
  PATH_SUFFIXES src/libscotch lib
  )

find_library(PTSCOTCH_LIBRARY_ESMUMPS ptesmumps
  PATHS ${PTSCOTCH_DIR}
  PATH_SUFFIXES src/libscotch lib
  )

find_path(PTSCOTCH_INCLUDE_PATH ptscotch.h
  PATHS ${PTSCOTCH_DIR}
  PATH_SUFFIXES include scotch src/libscotch include/scotch
  )

#===============================================================================
mark_as_advanced(PTSCOTCH_LIBRARY)
mark_as_advanced(PTSCOTCH_LIBRARY_ERR)
mark_as_advanced(PTSCOTCH_LIBRARY_ESMUMPS)
mark_as_advanced(PTSCOTCH_INCLUDE_PATH)

set(PTSCOTCH_LIBRARIES_ALL ${PTSCOTCH_LIBRARY} ${PTSCOTCH_LIBRARY_ERR})

if(PTSCOTCH_LIBRARY_ESMUMPS)
  set(PTSCOTCH_LIBRARIES_ALL ${PTSCOTCH_LIBRARY_ESMUMPS} ${PTSCOTCH_LIBRARIES_ALL})
endif()

set(PTSCOTCH_LIBRARIES ${PTSCOTCH_LIBRARIES_ALL} CACHE INTERNAL "Libraries for PT-Scotch" FORCE)

#===============================================================================
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PTSCOTCH DEFAULT_MSG
  PTSCOTCH_LIBRARY PTSCOTCH_LIBRARY_ERR PTSCOTCH_INCLUDE_PATH)


if(PTSCOTCH_INCLUDE_PATH)
  file(STRINGS ${PTSCOTCH_INCLUDE_PATH}/scotch.h PTSCOTCH_INCLUDE_CONTENT)
  string(REGEX MATCH "_cplusplus" _match ${PTSCOTCH_INCLUDE_CONTENT})
  if(_match)
    add_definitions(-DAKANTU_PTSCOTCH_NO_EXTERN)
  endif()
endif()

#===============================================================================
if(NOT PTSCOTCH_FOUND)
  set(PTSCOTCH_DIR "" CACHE PATH "Location of PT-Scotch library.")
endif(NOT PTSCOTCH_FOUND)
