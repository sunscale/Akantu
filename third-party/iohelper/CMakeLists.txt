#===============================================================================
# @file   CMakeLists.txt
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
#
# @date creation: Thu Oct 11 2012
# @date last modification: Fri Jun 13 2014
#
# @brief  main configuration file
#
# @section LICENSE
#
# Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#===============================================================================
# CMake Project
#===============================================================================
cmake_minimum_required(VERSION 2.6)

project(IOHelper)
enable_language(CXX)

#===============================================================================
# Misc.
#===============================================================================
set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries.")

find_package (ZLIB)

if(ZLIB_FOUND)
  set (IOHELPER_EXTERNAL_LIBS ${IOHELPER_EXTERNAL_LIBS} ${ZLIB_LIBRARIES})
  include_directories(${ZLIB_INCLUDE_DIR})
endif()

#===============================================================================
# Version Number
#===============================================================================
# IOHelper version number.  An even minor number corresponds to releases.
set(IOHELPER_MAJOR_VERSION 1)
set(IOHELPER_MINOR_VERSION 1)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/CMakeVersionGenerator.cmake)

define_project_version()


#==============================================================================
# Library
#===============================================================================
add_subdirectory(src)

#==============================================================================
# Packaging
#===============================================================================
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/IOHelperCPack.cmake)

#===============================================================================
# Tests
#===============================================================================
option(IOHELPER_TESTS "Activate tests" OFF)
if(IOHELPER_TESTS)
  enable_testing()
  include(CTest)
  add_subdirectory(test)
endif()
