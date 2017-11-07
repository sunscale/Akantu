#===============================================================================
# @file   CMakeVersionGenerator.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Thu Dec 20 2012
# @date last modification: Fri Jun 13 2014
#
# @brief  Set of macros used by akantu to handle the package system
#
# @section LICENSE
#
# Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# IOHelper is free  software: you can redistribute it and/or  modify it under the
# terms  of the  GNU Lesser  General Public  License as  published by  the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# IOHelper is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
# details.
#
# You should  have received  a copy  of the GNU  Lesser General  Public License
# along with IOHelper. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================

if(__DEFINE_PROJECT_VERSION__)
  return()
endif()
set(__DEFINE_PROJECT_VERSION__ TRUE)

macro(define_project_version)
  string(TOUPPER ${PROJECT_NAME} _project)

  if(EXISTS ${PROJECT_SOURCE_DIR}/VERSION)
    file(STRINGS ${PROJECT_SOURCE_DIR}/VERSION ${_project}_VERSION)

    if("${${_project}_VERSION}" MATCHES "^([0-9]+)")
      string(REGEX REPLACE "^([0-9]+).*" "\\1" _ver_major "${${_project}_VERSION}")
      set(${_project}_MAJOR_VERSION ${_ver_major})

      if("${${_project}_VERSION}" MATCHES "^${_ver_major}\\.([0-9]+)")
	string(REGEX REPLACE "^${_ver_major}\\.([0-9]+).*" "\\1" _ver_minor "${${_project}_VERSION}")
	set(${_project}_MINOR_VERSION ${_ver_minor})
	if("${${_project}_VERSION}" MATCHES "^${_ver_major}\\.${_ver_minor}\\.([0-9a-zA-Z\\-]+)")
	  string(REGEX REPLACE "^${_ver_major}\\.${_ver_minor}\\.([0-9a-zA-Z\\-]+).*" "\\1" _ver_build "${${_project}_VERSION}")
	  set(${_project}_BUILD_VERSION ${_ver_build})
	endif()
      endif()
    endif()
  # else()
  #   find_package(Subversion)

  #   if(SUBVERSION_FOUND)
  #     subversion_wc_info(${PROJECT_SOURCE_DIR} ${_project} ERROR_QUIET)
  #     if(${${_project}_WC_FOUND})
  #       set(${_project}_BUILD_VERSION ${${_project}_WC_REVISION})
  #       set(${_project}_VERSION
  #         "${${_project}_MAJOR_VERSION}.${${_project}_MINOR_VERSION}.${${_project}_BUILD_VERSION}"
  #         )
  #     endif()
  #   endif()
  endif()

  if(NOT ${_project}_VERSION)
    set(${_project}_VERSION
      "${${_project}_MAJOR_VERSION}.${${_project}_MINOR_VERSION}"
      )
  endif()

  # Append the library version information to the library target properties
  if(NOT ${_project}_NO_LIBRARY_VERSION)
    message(STATUS "${PROJECT_NAME} version: ${${_project}_VERSION}")

    set(${_project}_LIBRARY_PROPERTIES ${${_project}_LIBRARY_PROPERTIES}
      VERSION "${${_project}_VERSION}"
      SOVERSION "${${_project}_MAJOR_VERSION}.${${_project}_MINOR_VERSION}"
      )
  endif()
endmacro()