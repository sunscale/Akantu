#===============================================================================
# @file   iohelper.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Mar 30 2015
# @date last modification: Tue Jan 19 2016
#
# @brief  build script for iohelper
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

if(EXISTS ${PROJECT_SOURCE_DIR}/third-party/iohelper)
  set(IOHELPER_TARGETS_EXPORT ${AKANTU_TARGETS_EXPORT})
  add_subdirectory(${PROJECT_SOURCE_DIR}/third-party/iohelper)

  set(IOHELPER_SOURCE_TYPE "internal"
    CACHE INTERNAL "internal variable for clean export")

  #set(IOHELPER_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/third-party/iohelper/src" CACHE INTERNAL "IOHelper include directory")

  set(IOHELPER_LIBRARIES iohelper CACHE INTERNAL "IOHelper library")
  package_add_to_export_list(iohelper iohelper)
  mask_package_options(IOHELPER)
  return()
endif()



if(NOT EXISTS ${PROJECT_SOURCE_DIR}/third-party/${IOHELPER_ARCHIVE})
  set(_iohelper_download_command
    GIT_REPOSITORY ${IOHELPER_GIT}
    GIT_TAG ${IOHELPER_VERSION}
    )
else()
  set(_iohelper_download_command
    URL ${PROJECT_SOURCE_DIR}/third-party/${IOHELPER_ARCHIVE}
    )
endif()

if(CMAKE_VERSION VERSION_GREATER 3.1)
  set(_extra_options 
    UPDATE_DISCONNECTED 1
    DOWNLOAD_NO_PROGRESS 1
    EXCLUDE_FROM_ALL 1
    )
endif()

ExternalProject_Add(iohelper
  PREFIX ${PROJECT_BINARY_DIR}/third-party
  ${_iohelper_download_command}
  ${_extra_options}
  CMAKE_ARGS <SOURCE_DIR>/
  CMAKE_CACHE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE} -DCMAKE_CXX_COMPILER:PATH=${CMAKE_CXX_COMPILER}
  )

find_library(IOHELPER_LIBRARIES iohelper HINTS ${PROJECT_BINARY_DIR}/third-party)
if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
  set(_tmp ${IOHELPER_LIBRARIES})
  set(IOHELPER_LIBRARIES "${_tmp}.a" CACHE FILEPATH "" FORCE)
endif()

find_path(IOHELPER_INCLUDE_DIR iohelper/io_helper.hh HINTS ${PROJECT_BINARY_DIR}/third-party)

set(IOHELPER_SOURCE_TYPE "external" CACHE INTERNAL "internal variable for clean export")

mask_package_options(IOHELPER)

package_add_extra_dependency(IOHelper iohelper)

install(FILES ${IOHELPER_LIBRARIES}
  DESTINATION lib
  COMPONENT lib
  )

install(DIRECTORY ${IOHELPER_INCLUDE_DIR}
  DESTINATION include
  COMPONENT dev)
