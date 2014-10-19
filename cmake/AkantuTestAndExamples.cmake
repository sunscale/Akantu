#===============================================================================
# @file   AkantuTestAndExamples.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Oct 25 2010
# @date last modification: Tue Jun 24 2014
#
# @brief  macros for tests and examples
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
#
#===============================================================================

set(AKANTU_DIFF_SCRIPT ${AKANTU_CMAKE_DIR}/akantu_diff.sh)

#===============================================================================
function(manage_test_and_example et_name desc build_all label)
  string(TOUPPER ${et_name} upper_name)

  cmake_parse_arguments(manage_test_and_example
    ""
    "PACKAGE"
    ""
    ${ARGN}
    )

  set(_activated ON)
  if(manage_test_and_example_PACKAGE)
    list(FIND ${_project}_PACKAGE_SYSTEM_PACKAGES_ON ${manage_test_and_example_PACKAGE} _ret)
    if(_ret EQUAL -1)
      set(_activated OFF)
      file(RELATIVE_PATH _dir ${PROJECT_SOURCE_DIR}  ${CMAKE_CURRENT_SOURCE_DIR}/${et_name})
      list(APPEND AKANTU_TESTS_EXCLUDE_FILES /${_dir})
      set(AKANTU_TESTS_EXCLUDE_FILES ${AKANTU_TESTS_EXCLUDE_FILES} CACHE INTERNAL "")
    endif()
  endif()

  if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${et_name} AND _activated)
#    message("Example or test ${et_name} not present")
    return()
  endif()

  option(AKANTU_BUILD${label}${upper_name} "${desc}")
  mark_as_advanced(AKANTU_BUILD_${upper_name})

  if(${build_all} OR NOT _activated)
    set(AKANTU_BUILD${label}${upper_name}_OLD
      ${AKANTU_BUILD${label}${upper_name}}
      CACHE INTERNAL "${desc}" FORCE)

    set(AKANTU_BUILD${label}${upper_name} ${_activated}
      CACHE INTERNAL "${desc}" FORCE)
  else()
    if(DEFINED AKANTU_BUILD${label}${upper_name}_OLD)
      set(AKANTU_BUILD${label}${upper_name}
	${AKANTU_BUILD${label}${upper_name}_OLD}
	CACHE BOOL "${desc}" FORCE)

      unset(AKANTU_BUILD${label}${upper_name}_OLD
	CACHE)
    endif(DEFINED AKANTU_BUILD${label}${upper_name}_OLD)
  endif()

  if(AKANTU_BUILD${label}${upper_name})
    add_subdirectory(${et_name})
  endif(AKANTU_BUILD${label}${upper_name})
endfunction()

#===============================================================================
# Tests
#===============================================================================
if(AKANTU_TESTS)
  option(AKANTU_BUILD_ALL_TESTS "Build all tests")
#  mark_as_advanced(AKANTU_BUILD_ALL_TESTS)
endif(AKANTU_TESTS)

#===============================================================================
macro(register_test_old test_name)
  add_executable(${test_name} ${ARGN})
  target_link_libraries(${test_name} akantu ${AKANTU_EXTERNAL_LIBRARIES})

  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.sh)
    file(COPY ${test_name}.sh DESTINATION .)
  endif()

  if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${test_name}.sh)
    add_test(${test_name} ${CMAKE_CURRENT_BINARY_DIR}/${test_name}.sh)
  elseif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.verified)
    add_test(${test_name} ${AKANTU_DIFF_SCRIPT} ${test_name} ${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.verified)
  else()
    add_test(${test_name} ${CMAKE_CURRENT_BINARY_DIR}/${test_name})
  endif()
endmacro()

#===============================================================================
macro(add_akantu_test test_name desc)
  manage_test_and_example(${test_name} ${desc} AKANTU_BUILD_ALL_TESTS _ ${ARGN})
endmacro()


#===============================================================================
# Examples
#===============================================================================
if(AKANTU_EXAMPLES)
  option(AKANTU_BUILD_ALL_EXAMPLES "Build all examples")
#  mark_as_advanced(AKANTU_BUILD_ALL_EXAMPLES)
endif(AKANTU_EXAMPLES)

#===============================================================================
macro(register_example example_name)
  add_executable(${example_name} ${ARGN})
  target_link_libraries(${example_name} akantu ${AKANTU_EXTERNAL_LIBRARIES})
endmacro()

#===============================================================================
macro(add_example example_name desc)
  manage_test_and_example(${example_name} ${desc} AKANTU_BUILD_ALL_EXAMPLES _EXAMPLE_ ${ARGN})
endmacro()

#===============================================================================
function(register_test test_name)
  set(multi_variables
    SOURCES FILES_TO_COPY DEPENDENCIES DIRECTORIES_TO_CREATE COMPILE_OPTIONS EXTRA_FILES
    )

  cmake_parse_arguments(register_test
    ""
    "PACKAGE"
    "${multi_variables}"
    ${ARGN}
    )

  # add the test in a package if needed
  set(_test_define_package 0)
  set(_package_name AKANTU_CORE)
  if(register_test_PACKAGE)
    package_pkg_name(${register_test_PACKAGE} _package_name)
    set(_test_define_package 1)
    list(FIND ${_package_name}_TESTS ${test_name} _ret)
    if(_ret EQUAL -1)
      list(APPEND ${_package_name}_TESTS ${test_name})
    endif()
  endif()

  # check if the test should be activated
  set(_activate_test 0)
  foreach(_pkg ${${_project}_PACKAGE_SYSTEM_PACKAGES_ON})
    package_pkg_name(${_pkg} _package_name)
    list(FIND ${_package_name}_TESTS ${test_name} _ret)
    if(NOT _ret EQUAL -1)
      set(_activate_test 1)
    endif()
  endforeach()

  # check if the package is registered in at least a package
  set(_present_in_packages 0)
  foreach(_pkg ${${_project}_PACKAGE_SYSTEM_PACKAGES_NAMES_LIST_ALL})
    package_pkg_name(${_pkg} _package_name)
    list(FIND ${_package_name}_TESTS ${test_name} _ret)
    if(NOT _ret EQUAL -1)
      set(_present_in_packages 1)
    endif()
  endforeach()

  if(NOT _present_in_packages AND NOT _test_define_package)
    message("The test ${test_name} is not registered in any packages")
  endif()

  if(_activate_test)
    set(_source_file)
    foreach(_file ${register_test_SOURCES} ${register_test_UNPARSED_ARGUMENTS} ${register_test_EXTRA_FILES})
      if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${_file})
        list(APPEND _source_file ${CMAKE_CURRENT_SOURCE_DIR}/${_file})
      else()
        message("The file \"${_file}\" registred by the test \"${test_name}\" does not exists")
      endif()
    endforeach()

    add_executable(${test_name} ${register_test_SOURCES} ${register_test_UNPARSED_ARGUMENTS})
    set_property(TARGET ${test_name}  APPEND
      PROPERTY INCLUDE_DIRECTORIES ${AKANTU_INCLUDE_DIRS} ${AKANTU_EXTERNAL_LIB_INCLUDE_DIR})
    target_link_libraries(${test_name} akantu ${AKANTU_EXTERNAL_LIBRARIES})

    if(register_test_COMPILE_OPTIONS)
      set_target_properties(${test_name}
	PROPERTIES COMPILE_DEFINITIONS "${register_test_COMPILE_OPTIONS}")
    endif()

    if(register_test_FILES_TO_COPY)
      foreach(_file ${register_test_FILES_TO_COPY})
	file(COPY ${_file} DESTINATION .)
	list(APPEND _source_file ${CMAKE_CURRENT_SOURCE_DIR}/${_file})
      endforeach()
    endif()

    if(register_test_DIRECTORIES_TO_CREATE)
      foreach(_dir ${register_test_DIRECTORIES_TO_CREATE})
	if(IS_ABSOLUTE ${dir})
	  file(MAKE_DIRECTORY ${_dir})
	else()
  	  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${_dir})
	endif()
      endforeach()
    endif()

    foreach(_dep ${register_test_DEPENDENCIES})
      add_dependencies(${test_name} ${_dep})
      get_target_property(_dep_in_ressources ${_dep} RESSOURCES)
      if(_dep_in_ressources)
	list(APPEND _source_file ${_dep_in_ressources})
      endif()
    endforeach()

    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.sh)
      file(COPY ${test_name}.sh DESTINATION .)
      list(APPEND _source_file ${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.sh)
      add_test(${test_name}_run ${CMAKE_CURRENT_BINARY_DIR}/${test_name}.sh)
    elseif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.verified)
      list(APPEND _source_file ${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.verified)
      add_test(${test_name}_run ${AKANTU_DIFF_SCRIPT} ${test_name} ${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.verified)
    else()
      add_test(${test_name}_run ${CMAKE_CURRENT_BINARY_DIR}/${test_name})
    endif()

    set_tests_properties(${test_name}_run PROPERTIES DEPENDS ${test_name})

    set(_tmp ${AKANTU_TESTS_FILES})
    foreach(_file ${_source_file})
      get_filename_component(_full ${_file} ABSOLUTE)
      file(RELATIVE_PATH __file ${PROJECT_SOURCE_DIR} ${_full})
      list(APPEND _tmp ${__file})
      list(APPEND _pkg_tmp ${__file})
    endforeach()
    set(AKANTU_TESTS_FILES ${_tmp} CACHE INTERNAL "")
    set(${_package_name}_TESTS_FILES ${_pkg_tmp} CACHE INTERNAL "")
   endif()
endfunction()