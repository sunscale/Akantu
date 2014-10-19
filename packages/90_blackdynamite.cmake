#===============================================================================
# @file   blackdynamite.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date   Tue Nov 29 15:16:35 2011
#
# @brief  package description for BlackDynamite support
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
option(AKANTU_USE_BLACKDYNAMITE "Use blackdynamite library" OFF)
mark_as_advanced(AKANTU_USE_BLACKDYNAMITE)
find_package(Subversion)

if(SUBVERSION_FOUND)
  if(AKANTU_USE_BLACKDYNAMITE)
    set(BLACKDYNAMITE_SOURCE_DIR ${PROJECT_SOURCE_DIR}/third-party/blackdynamite)

    if(EXISTS ${BLACKDYNAMITE_SOURCE_DIR})
      execute_process(
	COMMAND ${Subversion_SVN_EXECUTABLE} up ${BLACKDYNAMITE_SOURCE_DIR}
	OUTPUT_VARIABLE _revision)
      string(REGEX REPLACE ".*At revision ([0-9]*)\\..*" "\\1" _rev "${_revision}")
      message(STATUS "Updating BlackDynamite: r${_rev} dynamite! dynamite!")
    else()
      message(STATUS "Checking out BlackDynamite: Can you digg it!")
      execute_process(
	COMMAND ${Subversion_SVN_EXECUTABLE} co svn+ssh://lsmssrv1.epfl.ch/space/repositories/SimulPack/BlackDynamite ${BLACKDYNAMITE_SOURCE_DIR}
	OUTPUT_QUIET)
    endif()


    set(BLACKDYNAMITE_TARGETS_EXPORT ${AKANTU_TARGETS_EXPORT})
    add_subdirectory(third-party/blackdynamite)

    list(APPEND AKANTU_EXTERNAL_LIBRARIES blackdynamite)
    list(APPEND AKANTU_EXTERNAL_LIB_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/third-party/blackdynamite/src ${BLACKDYNAMITE_EXTERNAL_INCLUDE_DIR})

    list(APPEND AKANTU_EXPORT_LIST blackdynamite)
    list(APPEND AKANTU_OPTION_LIST BLACKDYNAMITE)

    set(AKANTU_BLACKDYNAMITE_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/third-party/blackdynamite/src ${BLACKDYNAMITE_EXTERNAL_INCLUDE_DIR})

    set(AKANTU_BLACKDYNAMITE ON)
  else()
    set(AKANTU_BLACKDYNAMITE OFF)
  endif()
endif()
