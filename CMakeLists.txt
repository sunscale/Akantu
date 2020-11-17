#===============================================================================
# @file   CMakeLists.txt
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Jun 14 2010
# @date last modification: Fri Jan 22 2016
#
# @brief  main configuration file
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
# @section DESCRIPTION
#-------------------------------------------------------------------------------
#                      _               _
#                     | |             | |
#                 __ _| | ____ _ _ __ | |_ _   _
#                / _` | |/ / _` | '_ \| __| | | |
#               | (_| |   < (_| | | | | |_| |_| |
#                \__,_|_|\_\__,_|_| |_|\__|\__,_|
#
#===============================================================================

#===============================================================================
# CMake Project
#===============================================================================
cmake_minimum_required(VERSION 3.5.1)

# add this options before PROJECT keyword
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

if(CMAKE_VERSION VERSION_GREATER 3.12)
  cmake_policy(SET CMP0074 NEW)
endif()

set(AKANTU_MAJOR_VERSION 3)
set(AKANTU_MINOR_VERSION 1)
set(AKANTU_PATCH_VERSION 0)
set(AKANTU_COPYRIGHT "2010-2020, EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)")
set(AKANTU_MAINTAINER "Nicolas Richart")
set(AKANTU_HOMEPAGE_URL "htpps://akantu.ch")

if(CMAKE_VERSION VERSION_GREATER 3.12)
  project(Akantu
    HOMEPAGE_URL "htpps://akantu.ch")
else()
  project(Akantu)
endif()
  
enable_language(CXX)

#===============================================================================
# Misc. config for cmake
#===============================================================================
set(AKANTU_CMAKE_DIR "${PROJECT_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/Modules")

set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries.")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL
  "Enable/Disable output of compile commands during generation" FORCE)

mark_as_advanced(BUILD_SHARED_LIBS)

if(NOT AKANTU_TARGETS_EXPORT)
  set(AKANTU_TARGETS_EXPORT AkantuTargets)
endif()

include(CMakeVersionGenerator)
include(CMakePackagesSystem)
include(CMakeFlagsHandling)

include(AkantuPackagesSystem)
include(AkantuMacros)
include(AkantuCleaning)

#cmake_activate_debug_message()

include(GNUInstallDirs)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

#===============================================================================
# Version Number
#===============================================================================
# AKANTU version number.  An even minor number corresponds to releases.
define_project_version()

#===============================================================================
# Options
#===============================================================================
# Debug
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -DAKANTU_NDEBUG"
  CACHE STRING "Flags used by the compiler during release builds" FORCE)

#add_flags(cxx "-Wall -Wextra -pedantic -Werror")
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  add_flags(cxx "-Wall -Wextra -pedantic") # -Weffc++
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG_INIT} -ggdb3"
    CACHE STRING "Flags used by the compiler during debug builds" FORCE)
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT} -ggdb3"
    CACHE STRING "Flags used by the compiler during debug builds" FORCE)

  option (FORCE_COLORED_OUTPUT "Always produce ANSI-colored output (GNU/Clang only)." FALSE)
  mark_as_advanced(FORCE_COLORED_OUTPUT)
  if(FORCE_COLORED_OUTPUT)
    if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
      add_flags(cxx "-fcolor-diagnostics")
    else()
      add_flags(cxx "-fdiagnostics-color=always")
    endif()
  endif()
else()
  add_flags(cxx "-Wall")
endif()

option(AKANTU_EXAMPLES "Activate examples" OFF)
option(AKANTU_TESTS "Activate tests" OFF)

option(AKANTU_RUN_IN_DOCKER "Set the approriate flage tu run in docker" OFF)
set(AKANTU_PREFERRED_PYTHON_VERSION 3 CACHE STRING "Preferred version for python related things")
mark_as_advanced(AKANTU_PREFERRED_PYTHON_VERSION AKANTU_RUN_IN_DOCKER)

include(AkantuExtraCompilationProfiles)
#===============================================================================
# Dependencies
#===============================================================================
declare_akantu_types()

package_list_packages(${PROJECT_SOURCE_DIR}/packages
  EXTRA_PACKAGES_FOLDER ${PROJECT_SOURCE_DIR}/extra_packages
  NO_AUTO_COMPILE_FLAGS)

## meta option \todo better way to do it when multiple package give enable the
## same feature
if(AKANTU_SCOTCH)
  set(AKANTU_PARTITIONER ON)
else()
  set(AKANTU_PARTITIONER OFF)
endif()

if(AKANTU_MUMPS)
  set(AKANTU_SOLVER ON)
else()
  set(AKANTU_SOLVER OFF)
endif()

#===============================================================================
# Akantu library
#===============================================================================
add_subdirectory(src)

#===============================================================================
# Documentation
#===============================================================================
if(AKANTU_DOCUMENTATION_DEVELOPER_MANUAL OR AKANTU_DOCUMENTATION_USER_MANUAL)
  add_subdirectory(doc)
else()
  set(AKANTU_DOC_EXCLUDE_FILES "${PROJECT_SOURCE_DIR}/doc/manual" CACHE INTERNAL "")
endif()

#===============================================================================
# Python interface
#===============================================================================
package_is_activated(python_interface _python_act)
if(_python_act)
  if(IS_ABSOLUTE "${CMAKE_INSTALL_PREFIX}")
    set(AKANTU_PYTHON_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
  else()
    set(AKANTU_PYTHON_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_PREFIX}")
  endif()
  add_subdirectory(python)
endif()

#===============================================================================
# Examples and tests
#===============================================================================
include(AkantuTestsMacros)
include(AkantuExampleMacros)

if(AKANTU_TESTS)
  option(AKANTU_BUILD_ALL_TESTS "Build all tests" ON)
  find_package(GMSH REQUIRED)
endif()

# tests
add_test_tree(test)

if(AKANTU_EXAMPLES)
  if(AKANTU_TESTS)
    option(AKANTU_TEST_EXAMPLES "Run the examples" ON)
  endif()
  find_package(GMSH REQUIRED)
  add_subdirectory(examples)
endif()

#===============================================================================
# Install and Packaging
#===============================================================================
include(AkantuInstall)

option(AKANTU_DISABLE_CPACK
  "This option commands the generation of extra info for the \"make package\" target" ON)
mark_as_advanced(AKANTU_DISABLE_CPACK)
if(NOT AKANTU_DISABLE_CPACK)
  include(AkantuCPack)
endif()
