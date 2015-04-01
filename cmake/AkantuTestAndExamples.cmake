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
  option(AKANTU_BUILD_ALL_TESTS "Build all tests" ON)
  mark_as_advanced(AKANTU_BUILD_ALL_TESTS)
endif(AKANTU_TESTS)

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


# ==============================================================================
# this should be a macro due to the enable_testing
macro(add_test_tree dir)
  if(AKANTU_TESTS)
    enable_testing()
    include(CTest)
    mark_as_advanced(BUILD_TESTING)

    set(AKANTU_TESTS_EXCLUDE_FILES "" CACHE INTERNAL "")

    set(_akantu_current_parent_test ${dir} CACHE INTERNAL "Current test folder" FORCE)
    set(_akantu_${dir}_tests_count 0 CACHE INTERNAL "" FORCE)

    string(TOUPPER ${dir} _u_dir)
    set(AKANTU_BUILD_${_u_dir} ON CACHE INTERNAL "${desc}" FORCE)

    package_get_all_test_folders(_test_dirs)

    foreach(_dir ${_test_dirs})
      add_subdirectory(${_dir})
    endforeach()
  else()
    set(AKANTU_TESTS_EXCLUDE_FILES "${CMAKE_CURRENT_BINARY_DIR}/${dir}" CACHE INTERNAL "")
  endif()
endmacro()

# ==============================================================================
function(add_akantu_test dir desc)
  set(_my_parent_dir ${_akantu_current_parent_test})

  # initialize variables
  set(_akantu_current_parent_test ${dir} CACHE INTERNAL "Current test folder" FORCE)
  set(_akantu_${dir}_tests_count 0 CACHE INTERNAL "" FORCE)

  # set the option for this directory
  string(TOUPPER ${dir} _u_dir)
  option(AKANTU_BUILD_${_u_dir} "${desc}")
  mark_as_advanced(AKANTU_BUILD_${_u_dir})

  # add the sub-directory
  add_subdirectory(${dir})

  # if no test can be activated make the option disappear
  set(_force_deactivate_count FALSE)
  if(${_akantu_${dir}_tests_count} EQUAL 0)
    set(_force_deactivate_count TRUE)
  endif()

  # if parent off make the option disappear
  set(_force_deactivate_parent FALSE)
  string(TOUPPER ${_my_parent_dir} _u_parent_dir)
  if(NOT AKANTU_BUILD_${_u_parent_dir})
    set(_force_deactivate_parent TRUE)
  endif()

  if(_force_deactivate_parent OR _force_deactivate_count OR AKANTU_BUILD_ALL_TESTS)
    if(NOT DEFINED _AKANTU_BUILD_${_u_dir}_SAVE)
      set(_AKANTU_BUILD_${_u_dir}_SAVE ${AKANTU_BUILD_${_u_dir}} CACHE INTERNAL "" FORCE)
    endif()
    unset(AKANTU_BUILD_${_u_dir} CACHE)
    if(AKANTU_BUILD_ALL_TESTS AND NOT _force_deactivate_count)
      set(AKANTU_BUILD_${_u_dir} ON CACHE INTERNAL "${desc}" FORCE)
    else()
      set(AKANTU_BUILD_${_u_dir} OFF CACHE INTERNAL "${desc}" FORCE)
    endif()
  else()
    if(DEFINED _AKANTU_BUILD_${_u_dir}_SAVE)
      unset(AKANTU_BUILD_${_u_dir} CACHE)
      set(AKANTU_BUILD_${_u_dir} ${_AKANTU_BUILD_${_u_dir}_SAVE} CACHE BOOL "${desc}")
      unset(_AKANTU_BUILD_${_u_dir}_SAVE CACHE)
    endif()
  endif()

  # adding up to the parent count
  math(EXPR _tmp_parent_count "${_akantu_${dir}_tests_count} + ${_akantu_${_my_parent_dir}_tests_count}")
  set(_akantu_${_my_parent_dir}_tests_count ${_tmp_parent_count} CACHE INTERNAL "" FORCE)

  # restoring the parent current dir
  set(_akantu_current_parent_test ${_my_parent_dir} CACHE INTERNAL "Current test folder" FORCE)
endfunction()


# ==============================================================================
function(register_test test_name)
   set(multi_variables
    SOURCES FILES_TO_COPY DEPENDENCIES DIRECTORIES_TO_CREATE COMPILE_OPTIONS EXTRA_FILES
    )

  cmake_parse_arguments(_register_test
    ""
    "PACKAGE"
    "${multi_variables}"
    ${ARGN}
    )

  if(NOT _register_test_PACKAGE)
    message(FATAL_ERROR "No reference package was defined for the test ${test_name} in folder ${CMAKE_CURRENT_SOURCE_DIR}")
  endif()

  package_is_activated(${_register_test_PACKAGE} _act)

  if(_act)
    math(EXPR _tmp_parent_count "${_akantu_${_akantu_current_parent_test}_tests_count} + 1")
    set(_akantu_${_akantu_current_parent_test}_tests_count ${_tmp_parent_count} CACHE INTERNAL "" FORCE)

    string(TOUPPER ${_akantu_current_parent_test} _u_parent)
    if(AKANTU_BUILD_${_u_parent})
      package_get_all_include_directories(
	AKANTU_LIBRARY_INCLUDE_DIRS
	)

      package_get_all_external_informations(
	AKANTU_EXTERNAL_INCLUDE_DIR
	AKANTU_EXTERNAL_LIBRARIES
	)
      # set the proper includes to build most of the tests
      include_directories(
	${AKANTU_INCLUDE_DIRS}
	${AKANTU_EXTERNAL_LIB_INCLUDE_DIR}
	)

      # Register the executable to compile
      add_executable(${test_name} ${_register_test_SOURCES} ${_register_test_UNPARSED_ARGUMENTS})
      set_property(TARGET ${test_name}  APPEND
	PROPERTY INCLUDE_DIRECTORIES ${AKANTU_LIBRARY_INCLUDE_DIRS} ${AKANTU_EXTERNAL_INCLUDE_DIR})
      target_link_libraries(${test_name} akantu ${AKANTU_EXTERNAL_LIBRARIES})

      # add the extra compilation options
      if(_register_test_COMPILE_OPTIONS)
	set_target_properties(${test_name}
	  PROPERTIES COMPILE_DEFINITIONS "${_register_test_COMPILE_OPTIONS}")
      endif()

      set(_test_all_files)
      # add the different dependencies (meshes, local libraries, ...)
      foreach(_dep ${_register_test_DEPENDENCIES})
	add_dependencies(${test_name} ${_dep})
	get_target_property(_dep_in_ressources ${_dep} RESSOURCES)

	if(_dep_in_ressources)
	  list(APPEND _test_all_files ${_dep_in_ressources})
	endif()
      endforeach()

      # copy the needed files to the build folder
      if(_register_test_FILES_TO_COPY)
	foreach(_file ${_register_test_FILES_TO_COPY})
	  file(COPY ${_file} DESTINATION .)
	  list(APPEND _test_all_files ${CMAKE_CURRENT_SOURCE_DIR}/${_file})
	endforeach()
      endif()

      # create the needed folders in the build folder
      if(_register_test_DIRECTORIES_TO_CREATE)
	foreach(_dir ${_register_test_DIRECTORIES_TO_CREATE})
	  if(IS_ABSOLUTE ${dir})
	    file(MAKE_DIRECTORY ${_dir})
	  else()
  	    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${_dir})
	  endif()
	endforeach()
      endif()

      # add the source files in the list of all files
      foreach(_file ${_register_test_SOURCES} ${_register_test_UNPARSED_ARGUMENTS} ${_register_test_EXTRA_FILES})
	if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${_file})
          list(APPEND _test_all_files ${CMAKE_CURRENT_SOURCE_DIR}/${_file})
	else()
          message("The file \"${_file}\" registred by the test \"${test_name}\" does not exists")
	endif()
      endforeach()

      # register the test for ctest
      if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.sh)
	file(COPY ${test_name}.sh DESTINATION .)
	list(APPEND _test_all_files ${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.sh)
	add_test(${test_name}_run
	  ${CMAKE_CURRENT_BINARY_DIR}/${test_name}.sh)
      elseif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.verified)
	list(APPEND _test_all_files ${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.verified)
	add_test(${test_name}_run
	  ${AKANTU_DIFF_SCRIPT} ${test_name} ${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.verified)
      else()
	add_test(${test_name}_run
	  ${CMAKE_CURRENT_BINARY_DIR}/${test_name})
      endif()

      # add the executable as a dependency of the run
      set_tests_properties(${test_name}_run PROPERTIES DEPENDS ${test_name})

      # clean the list of all files for this test and add them in the total list
      set(_tmp ${AKANTU_TESTS_FILES})
      foreach(_file ${_source_file})
	get_filename_component(_full ${_file} ABSOLUTE)
	file(RELATIVE_PATH __file ${PROJECT_SOURCE_DIR} ${_full})
	list(APPEND _tmp ${__file})
	list(APPEND _pkg_tmp ${__file})
      endforeach()
      set(AKANTU_TESTS_FILES ${_tmp} CACHE INTERNAL "")
    endif()
  endif()
endfunction()