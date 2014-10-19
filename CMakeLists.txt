#===============================================================================
# @file   CMakeLists.txt
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Jun 14 2010
# @date last modification: Mon Sep 15 2014
#
# @brief  main configuration file
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
cmake_minimum_required(VERSION 2.8)

# add this options before PROJECT keyword
#set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project(Akantu)
enable_language(CXX)

#===============================================================================
# Misc. config for cmake
#===============================================================================
set(AKANTU_CMAKE_DIR "${PROJECT_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/Modules")
set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries.")
mark_as_advanced(BUILD_SHARED_LIBS)

if(NOT AKANTU_TARGETS_EXPORT)
  set(AKANTU_TARGETS_EXPORT AkantuLibraryDepends)
endif()

include(CMakeVersionGenerator)
include(CMakePackagesSystem)
include(CMakeFlagsHandling)

include(AkantuMacros)

#cmake_activate_debug_message()

#===============================================================================
# Version Number
#===============================================================================
# AKANTU version number.  An even minor number corresponds to releases.
set(AKANTU_MAJOR_VERSION 2)
set(AKANTU_MINOR_VERSION 2)
set(AKANTU_PATCH_VERSION 0)
define_project_version()

#===============================================================================
# Options
#===============================================================================
# Debug
#option(AKANTU_DEBUG "Compiles akantu with the debug messages" ON)
#mark_as_advanced(AKANTU_DEBUG)

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -DAKANTU_NDEBUG"
  CACHE STRING "Flags used by the compiler during release builds" FORCE)

add_flags(cxx "-Wall")

# if(NOT AKANTU_DEBUG)
#   set(AKANTU_NDEBUG ON)
# else()
#   set(AKANTU_NDEBUG OFF)
# endif()

#Profiling
set(CMAKE_CXX_FLAGS_PROFILING "-g -pg" CACHE STRING "Flags used by the compiler during profiling builds")
set(CMAKE_EXE_LINKER_FLAGS_PROFILING "-pg" CACHE STRING "Flags used by the linker during profiling builds")
set(CMAKE_SHARED_LINKER_FLAGS_PROFILING "-pg" CACHE STRING "Flags used by the linker during profiling builds")
mark_as_advanced(CMAKE_CXX_FLAGS_PROFILING)
mark_as_advanced(CMAKE_EXE_LINKER_FLAGS_PROFILING)
mark_as_advanced(CMAKE_SHARED_LINKER_FLAGS_PROFILING)

include(CMakeDetermineCCompiler)
cmake_deactivate_debug_message(PackagesSystem)
#===============================================================================
# Dependencies
#===============================================================================
set(AKANTU_BOOST_COMPONENT)
set(AKANTU_EXTRA_TARGET_DEPENDENCIES)
add_all_packages(${PROJECT_SOURCE_DIR}/packages ${PROJECT_SOURCE_DIR}/src)

include_boost()

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
if(AKANTU_DOCUMENTATION_DOXYGEN OR AKANTU_DOCUMENTATION_MANUAL)
  add_subdirectory(doc)
else()
  set(AKANTU_DOC_EXCLUDE_FILES "${PROJECT_SOURCE_DIR}/doc/manual" CACHE INTERNAL "")
endif()

#===============================================================================
# Examples and tests
#===============================================================================
option(AKANTU_EXAMPLES "Activate examples" OFF)
option(AKANTU_TESTS "Activate tests" OFF)

if(AKANTU_EXAMPLES OR AKANTU_TESTS)
  include(AkantuTestAndExamples)
  find_package(GMSH REQUIRED)
endif()

if(AKANTU_EXAMPLES)
  add_subdirectory(examples)
endif()

if(AKANTU_TESTS)
  enable_testing()
  include(CTest)
  mark_as_advanced(BUILD_TESTING)
  set(AKANTU_TESTS_EXCLUDE_FILES "" CACHE INTERNAL "")
  add_subdirectory(test)
else()
  set(AKANTU_TESTS_EXCLUDE_FILES "${PROJECT_SOURCE_DIR}/test" CACHE INTERNAL "")
endif()


#===============================================================================
# Install and Packaging
#===============================================================================
include(AkantuInstall)

if(NOT AKANTU_DISABLE_CPACK)
  include(AkantuCPack)
endif()
