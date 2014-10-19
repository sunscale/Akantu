#===============================================================================
# @file   cpp_array.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date   Mon Nov 21 18:19:15 2011
#
# @brief  package description for cpp_array project
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
option(AKANTU_USE_CPPARRAY "Use cpp-array library" OFF)
option(AKANTU_USE_THIRD_PARTY_CPPARRAY "Automatic download of the CPP-ARRAY library" ON)
mark_as_advanced(AKANTU_USE_THIRD_PARTY_CPPARRAY AKANTU_USE_CPPARRAY)

if(AKANTU_USE_CPPARRAY AND AKANTU_USE_THIRD_PARTY_CPPARRAY)
  find_package(Git)

  if(GIT_FOUND)
    if(EXISTS ${PROJECT_SOURCE_DIR}/third-party/cpp-array)
      execute_process(
	COMMAND ${GIT_EXECUTABLE} pull
	WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/third-party/cpp-array
	OUTPUT_VARIABLE _revision)
      message(STATUS "Updating Cpp-Array")
    else()
      message(STATUS "Cloning Cpp-Array")
      execute_process(
	COMMAND ${GIT_EXECUTABLE} clone https://code.google.com/p/cpp-array.git ${PROJECT_SOURCE_DIR}/third-party/cpp-array
	OUTPUT_QUIET)
    endif()
  endif()

  if(EXISTS ${PROJECT_SOURCE_DIR}/third-party/cpp-array/)
    set(cpp-array_TESTS OFF CACHE BOOL "cpparray tests" FORCE)
    add_subdirectory(${PROJECT_SOURCE_DIR}/third-party/cpp-array/)
    set(cpp-array_TESTS OFF CACHE BOOL "cpparray tests" FORCE)

    mark_as_advanced(cpp-array_DEV)
    mark_as_advanced(cpp-array_DOCUMENTATION)
    mark_as_advanced(cpp-array_TESTS)
    mark_as_advanced(CUDA)
    mark_as_advanced(ARRAY_USER_LIB_PATH)

    list(APPEND AKANTU_EXTERNAL_LIB_INCLUDE_DIR ${cpp-array_INCLUDE_DIRS} ${CPP-ARRAY_INCLUDE_DIRS})
    list(APPEND AKANTU_EXTERNAL_LIBRARIES ${CPP-ARRAY_LIBRARIES})
    list(APPEND AKANTU_EXTERNAL_LIB_INCLUDE_DIR ${CPP-ARRAY_INCLUDE_DIRS})
    list(APPEND CPACK_SOURCE_IGNORE_FILES ${PROJECT_SOURCE_DIR}/third-party/cpp-array/)
    set(AKANTU_CPPARRAY_INCLUDE_DIR ${cpp-array_INCLUDE_DIRS} ${CPP-ARRAY_INCLUDE_DIRS})


    set(AKANTU_CPPARRAY ON)
    list(APPEND AKANTU_OPTION_LIST CPPARRAY)
    set(AKANTU_CPPARRAY ${AKANTU_CPPARRAY} CACHE INTERNAL "Use cpp-array library" FORCE)
  else()
    message(STATUS "Cpp-Array could not be found! Please install git and/or place cpp-array in the third-party folder of Akantu")
  endif()
else()
  add_optional_external_package(CppArray "Use cpp-array library" OFF)
endif()


set(AKANTU_CPPARRAY_DOCUMENTATION "
This package provides access to the \\href{https://code.google.com/p/cpp-array/}{cpp-array} open-source project. If internet is accessible when configuring the project (during cmake call) this package will be auto-downloaded.
")
