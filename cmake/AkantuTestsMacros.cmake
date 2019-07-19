#===============================================================================
# @file   AkantuTestsMacros.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Sep 03 2010
# @date last modification: Fri Jan 22 2016
#
# @brief  macros for tests
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
#===============================================================================

#[=======================================================================[.rst:
AkantuTestsMacros
-----------------

This modules provides the functions to helper to declare tests and folders
containing tests in akantu

.. command:: add_test_tree

add_test_tree(<test_direcotry>)

``<test_directory>`` is the entry direcroty of the full structure of
subfolders containing tests

.. command:: add_akantu_test

add_akantu_test(<dir> <desc>)

This function add a subdirectory ``<dir>`` of tests that will be conditionnaly
activable and will be visible only if the parent folder as been activated An
option ``AKANTU_BUILD_TEST_<dir>`` will appear in ccmake with the description
``<desc>``. The compilation of all tests can be forced with the option
``AKANTU_BUILD_ALL_TESTS``

.. command:: register_test

register_test(<test_name>
  SOURCES <sources>...
  PACKAGE <akantu_packages>...
  SCRIPT <scirpt>
  [FILES_TO_COPY <filenames>...]
  [DEPENDS <targets>...]
  [DIRECTORIES_TO_CREATE <directories>...]
  [COMPILE_OPTIONS <flags>...]
  [EXTRA_FILES <filnames>...]
  [LINK_LIBRARIES <libraries>...]
  [INCLUDE_DIRECTORIES <include>...]
  [UNSABLE]
  [PARALLEL]
  [PARALLEL_LEVEL <procs>...]
  )

This function defines a test ``<test_name>_run`` this test could be of
different nature depending on the context. If Just sources are provided the
test consist of running the executable generated. If a file ``<test_name>.sh``
is present the test will execute the script. And if a ``<test_name>.verified``
exists the output of the test will be compared to this reference file

The options are:

``SOURCES <sources>...``
The list of source files to compile to generate the executable of the test

``PACKAGE <akantu_packages>...``
The list of package to which this test belongs. The test will be activable
only of all the packages listed are activated

``SCRIPT <script>``
The script to execute instead of the executable

``FILES_TO_COPY <filenames>...``
List of files to copy from the source directory to the build directory

``DEPENDS <targets>...``
List of targets the test depends on, for example if a mesh as to be generated

``DIRECTORIES_TO_CREATE <directories>...``
Obsolete. This specifies a list of directories that have to be created in
the build folder

``COMPILE_OPTIONS <flags>...``
List of extra compilations options to pass to the compiler

``EXTRA_FILES <filnames>...``
Files to consider when generating a package_source

``UNSABLE``
If this option is specified the test can be unacitivated by the glocal option
``AKANTU_BUILD_UNSTABLE_TESTS``, this is mainly intendeed to remove test
under developement from the continious integration

``PARALLEL``
This specifies that this test should be run in parallel. It will generate a
series of test for different number of processors. This automaticaly adds a
dependency to the package ``AKANTU_PARALLEL``

``PARALLEL_LEVEL``
This defines the different processor numbers to use, if not defined the
macro tries to determine it in a "clever" way

]=======================================================================]

set(AKANTU_DRIVER_SCRIPT ${AKANTU_CMAKE_DIR}/akantu_test_driver.sh)

# ==============================================================================
macro(add_test_tree dir)
  if(AKANTU_TESTS)
    enable_testing()
    include(CTest)
    mark_as_advanced(BUILD_TESTING)

    set(_akantu_current_parent_test ${dir} CACHE INTERNAL "Current test folder" FORCE)
    set(_akantu_${dir}_tests_count 0 CACHE INTERNAL "" FORCE)

    string(TOUPPER ${dir} _u_dir)
    set(AKANTU_BUILD_${_u_dir} ON CACHE INTERNAL "${desc}" FORCE)

    package_get_all_test_folders(_test_dirs)

    foreach(_dir ${_test_dirs})
      add_subdirectory(${_dir})
    endforeach()
  endif()
endmacro()


set(_test_flags
  UNSTABLE
  PARALLEL
  PYTHON
  GTEST
  HEADER_ONLY
  )

set(_test_one_variables
  POSTPROCESS
  SCRIPT
  )

set(_test_multi_variables
  SOURCES
  FILES_TO_COPY
  DEPENDS
  DIRECTORIES_TO_CREATE
  COMPILE_OPTIONS
  EXTRA_FILES
  LINK_LIBRARIES
  INCLUDE_DIRECTORIES
  PACKAGE
  PARALLEL_LEVEL
  )


# ==============================================================================
function(add_akantu_test dir desc)
  if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${dir})
    return()
  endif()

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

function(is_test_active is_active)
  cmake_parse_arguments(_register_test
    "${_test_flags}"
    "${_test_one_variables}"
    "${_test_multi_variables}"
    ${ARGN}
    )

  if(NOT _register_test_PACKAGE)
    message(FATAL_ERROR "No reference package was defined for the test"
      " ${test_name} in folder ${CMAKE_CURRENT_SOURCE_DIR}")
  endif()

  if(_register_test_PYTHON)
    list(APPEND _register_test_PACKAGE python_interface)
  endif()

  set(_test_act TRUE)
  # Activate the test anly if all packages associated to the test are activated
  foreach(_package ${_register_test_PACKAGE})
    package_is_activated(${_package} _act)
    if(NOT _act)
      set(_test_act FALSE)
    endif()
  endforeach()

  # check if the test is marked unstable and if the unstable test should be run
  if(_register_test_UNSTABLE AND NOT AKANTU_BUILD_UNSTABLE_TESTS)
    set(_test_act FALSE)
  endif()

  if(_test_act)
    # todo this should be checked for the build package_sources since the file will not be listed.
    math(EXPR _tmp_parent_count "${_akantu_${_akantu_current_parent_test}_tests_count} + 1")
    set(_akantu_${_akantu_current_parent_test}_tests_count ${_tmp_parent_count} CACHE INTERNAL "" FORCE)
  endif()

  string(TOUPPER ${_akantu_current_parent_test} _u_parent)
  if(NOT (AKANTU_BUILD_${_u_parent} OR AKANTU_BUILD_ALL_TESTS))
    set(_test_act FALSE)
  endif()

  set(${is_active} ${_test_act} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
function(register_gtest_sources)
  cmake_parse_arguments(_register_test
    "${_test_flags}"
    "${_test_one_variables}"
    "${_test_multi_variables}"
    ${ARGN}
    )

  is_test_active(_is_active ${ARGN})
  register_test_files_to_package(${ARGN})

  if(NOT _is_active)
    return()
  endif()

  if(_register_test_PACKAGE)
    set(_list ${_gtest_PACKAGE})
    list(APPEND _list ${_register_test_PACKAGE})
    list(REMOVE_DUPLICATES _list)
    set(_gtest_PACKAGE ${_list} PARENT_SCOPE)
  endif()

  foreach (_var ${_test_flags})
    if(_var STREQUAL "HEADER_ONLY")
      if(NOT DEFINED_register_test_${_var})
        set(_gtest_${_var} OFF PARENT_SCOPE)
      elseif(NOT DEFINED _gtest_${_var})
        set(_gtest_${_var} ON PARENT_SCOPE)
      endif()
      continue()
    endif()

    if(_register_test_${_var})
      set(_gtest_${_var} ON PARENT_SCOPE)
    else()
      if(_gtest_${_var})
        message("Another gtest file required ${_var} to be ON it will be globally set for this folder...")
      endif()
    endif()
  endforeach()

  if(_register_test_UNPARSED_ARGUMENTS)
    list(APPEND _register_test_SOURCES ${_register_test_UNPARSED_ARGUMENTS})
  endif()

  foreach (_var ${_test_multi_variables})
    if(_register_test_${_var})
      set(_list ${_gtest_${_var}})
      list(APPEND _list ${_register_test_${_var}})
      list(REMOVE_DUPLICATES _list)
      set(_gtest_${_var} ${_list} PARENT_SCOPE)
    endif()
  endforeach()
endfunction()

# ==============================================================================
function(akantu_pybind11_add_module target)
  package_is_activated(pybind11 _pybind11_act)
  if(_pybind11_act)
    package_get_all_external_informations(
      INTERFACE_INCLUDE AKANTU_INTERFACE_EXTERNAL_INCLUDE_DIR
      )

    pybind11_add_module(${target} ${ARGN})
    target_include_directories(${target} SYSTEM PRIVATE ${PYBIND11_INCLUDE_DIR}
      ${AKANTU_INTERFACE_EXTERNAL_INCLUDE_DIR})
    set_property(TARGET ${target} PROPERTY DEBUG_POSTFIX "")
  endif()
endfunction()

# ==============================================================================
function(register_gtest_test test_name)
  if(NOT _gtest_PACKAGE)
    return()
  endif()

  set(_argn ${test_name}_gtest)

  set(_link_libraries GTest::GTest GTest::Main)

  list(FIND _gtest_PACKAGE python_interface _pos)
  package_is_activated(python_interface _python_interface_act)

  if(_python_interface_act AND (NOT _pos EQUAL -1))
    list(APPEND _link_libraries pybind11::embed)
    set(_compile_flags COMPILE_OPTIONS "AKANTU_TEST_USE_PYBIND11")
  endif()

  is_test_active(_is_active ${ARGN} PACKAGE ${_gtest_PACKAGE})
  if(NOT _is_active)
    return()
  endif()

  register_gtest_sources(${ARGN}
    SOURCES ${PROJECT_SOURCE_DIR}/test/test_gtest_main.cc
    LINK_LIBRARIES ${_link_libraries}
    PACKAGE ${_gtest_PACKAGE}
    ${_compile_flags}
    )

  foreach (_var ${_test_flags})
    if(_gtest_${_var})
      list(APPEND _argn ${_var})
      unset(_gtest_${_var})
    endif()
  endforeach()

  foreach (_var ${_test_multi_variables})
    if(_gtest_${_var})
      list(APPEND _argn ${_var} ${_gtest_${_var}})
      unset(_gtest_${_var})
    endif()
  endforeach()

  register_test(${_argn} GTEST)
  target_include_directories(${test_name}_gtest PRIVATE ${PROJECT_SOURCE_DIR}/test)
endfunction()

# ==============================================================================
function(register_test test_name)
  cmake_parse_arguments(_register_test
    "${_test_flags}"
    "${_test_one_variables}"
    "${_test_multi_variables}"
    ${ARGN}
    )

  register_test_files_to_package(${ARGN})
  is_test_active(_test_act ${ARGN})
  if(NOT _test_act)
    return()
  endif()

  set(_extra_args)

  # check that the sources are files that need to be compiled
  if(_register_test_SOURCES} OR _register_test_UNPARSED_ARGUMENTS)
    set(_need_to_compile TRUE)
  else()
    set(_need_to_compile FALSE)
  endif()

  set(_compile_source)
  foreach(_file ${_register_test_SOURCES} ${_register_test_UNPARSED_ARGUMENTS})
    if(_file MATCHES "\\.cc$" OR _file MATCHES "\\.hh$")
      list(APPEND _compile_source ${_file})
    endif()
  endforeach()

  if(_compile_source)
    # get the include directories for sources in activated directories
    package_get_all_include_directories(
      AKANTU_LIBRARY_INCLUDE_DIRS
      )

    # get the external packages compilation and linking informations
    package_get_all_external_informations(
      INTERFACE_INCLUDE AKANTU_EXTERNAL_INCLUDE_DIR
      )

    foreach(_pkg ${_register_test_PACKAGE})
      package_get_nature(${_pkg} _nature)
      if(_nature MATCHES "^external.*")
        package_get_include_dir(${_pkg} _incl)
        package_get_libraries(${_pkg} _libs)

        list(APPEND _register_test_INCLUDE_DIRECTORIES ${_incl})
        list(APPEND _register_test_LINK_LIBRARIES ${_libs})
      endif()
    endforeach()

    # Register the executable to compile
    add_executable(${test_name} ${_compile_source})

    # set the proper includes to build most of the tests
    target_include_directories(${test_name}
      PRIVATE ${AKANTU_LIBRARY_INCLUDE_DIRS}
              ${AKANTU_EXTERNAL_INCLUDE_DIR}
              ${PROJECT_BINARY_DIR}/src
              ${_register_test_INCLUDE_DIRECTORIES})

    if(NOT _register_test_HEADER_ONLY)
      target_link_libraries(${test_name} PRIVATE akantu ${_register_test_LINK_LIBRARIES})
    else()
      get_target_property(_features akantu INTERFACE_COMPILE_FEATURES)
      target_link_libraries(${test_name} ${_register_test_LINK_LIBRARIES})
      target_compile_features(${test_name} PRIVATE ${_features})
    endif()

    # add the extra compilation options
    if(_register_test_COMPILE_OPTIONS)
      set_target_properties(${test_name}
        PROPERTIES COMPILE_DEFINITIONS "${_register_test_COMPILE_OPTIONS}")
    endif()

    if(AKANTU_EXTRA_CXX_FLAGS)
      set_target_properties(${test_name}
        PROPERTIES COMPILE_FLAGS "${AKANTU_EXTRA_CXX_FLAGS}")
    endif()
  else()
    add_custom_target(${test_name} ALL)
    if(_register_test_UNPARSED_ARGUMENTS AND NOT _register_test_SCRIPT)
      set(_register_test_SCRIPT ${_register_test_UNPARSED_ARGUMENTS})
    endif()
  endif()

  if(_register_test_DEPENDS)
    add_dependencies(${test_name} ${_register_test_DEPENDS})
  endif()

  # copy the needed files to the build folder
  if(_register_test_FILES_TO_COPY)
    foreach(_file ${_register_test_FILES_TO_COPY})
      _add_file_to_copy(${test_name} "${_file}")
    endforeach()
  endif()

  # create the needed folders in the build folder
  if(_register_test_DIRECTORIES_TO_CREATE)
    foreach(_dir ${_register_test_DIRECTORIES_TO_CREATE})
      if(IS_ABSOLUTE ${dir})
        file(MAKE_DIRECTORY "${_dir}")
      else()
        file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${_dir}")
      endif()
    endforeach()
  endif()

  # register the test for ctest
  set(_arguments -n "${test_name}")
  if(_register_test_SCRIPT)
    _add_file_to_copy(${test_name} ${_register_test_SCRIPT})
    if(_register_test_PYTHON)
      if(NOT PYTHONINTERP_FOUND)
        find_package(PythonInterp ${AKANTU_PREFERRED_PYTHON_VERSION} REQUIRED)
      endif()
      list(APPEND _arguments -e "${PYTHON_EXECUTABLE}")
      list(APPEND _extra_args "${_register_test_SCRIPT}")
      add_dependencies(${test_name} py11_akantu)
    else()
      list(APPEND _arguments -e "./${_register_test_SCRIPT}")
    endif()
  elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.sh")
    _add_file_to_copy(${test_name} ${test_name}.sh)
    list(APPEND _arguments -e "./${test_name}.sh")
  else()
    list(APPEND _arguments -e "./${test_name}")
  endif()

  if(_register_test_GTEST)
    list(APPEND _extra_args "--" "--gtest_output=xml:${PROJECT_BINARY_DIR}/gtest_reports/${test_name}.xml")
  endif()
  
  list(APPEND _arguments -E "${PROJECT_BINARY_DIR}/akantu_environement.sh")

  package_is_activated(parallel _is_parallel)
  if(_is_parallel AND AKANTU_TESTS_ALWAYS_USE_MPI AND NOT _register_test_PARALLEL)
    set(_register_test_PARALLEL TRUE)
    set(_register_test_PARALLEL_LEVEL 1)
  endif()

  if(_register_test_PARALLEL AND _is_parallel)
    set(_exe ${MPIEXEC})
    if(NOT _exe)
      set(_exe ${MPIEXEC_EXECUTABLE})
    endif()
    list(APPEND _arguments -p "${_exe} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG}")
    if(_register_test_PARALLEL_LEVEL)
      set(_procs "${_register_test_PARALLEL_LEVEL}")
    elseif(CMAKE_VERSION VERSION_GREATER "3.0")
      set(_procs)
      include(ProcessorCount)
      ProcessorCount(N)
      while(N GREATER 1)
        list(APPEND _procs ${N})
        math(EXPR N "${N} / 2")
      endwhile()
      list(APPEND _procs 1)
    endif()

    if(NOT _procs)
      set(_procs 2)
    endif()
  endif()

  if(_register_test_POSTPROCESS)
    list(APPEND _arguments -s "${_register_test_POSTPROCESS}")
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${_register_test_POSTPROCESS}
      FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  endif()

  list(APPEND _arguments -w "${CMAKE_CURRENT_BINARY_DIR}")

  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.verified")
    list(APPEND _arguments -r "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.verified")
  endif()

  string(REPLACE ";" " " _command "${_arguments}")

  # register them test
  if(_procs)
    foreach(p ${_procs})
      add_test(NAME ${test_name}_${p} COMMAND ${AKANTU_DRIVER_SCRIPT} ${_arguments} -N ${p} ${_extra_args})
      set_property(TEST ${test_name}_${p} PROPERTY PROCESSORS ${p})
    endforeach()
  else()
    add_test(NAME ${test_name} COMMAND ${AKANTU_DRIVER_SCRIPT} ${_arguments} ${_extra_args})
    set_property(TEST ${test_name} PROPERTY PROCESSORS 1)
  endif()
endfunction()


function(register_test_files_to_package)
  cmake_parse_arguments(_register_test
    "${_test_flags}"
    "${_test_one_variables}"
    "${_test_multi_variables}"
    ${ARGN}
    )

  if(_register_test_PYTHON)
    list(APPEND _register_test_PACKAGE python_interface)
  endif()

  set(_test_all_files)
  # add the source files in the list of all files
  foreach(_file ${_register_test_SOURCES} ${_register_test_UNPARSED_ARGUMENTS}
      ${_register_test_EXTRA_FILES} ${_register_test_SOURCES} ${_register_test_SCRIPT}
      ${_register_test_POSTPROCESS} ${_register_test_FILES_TO_COPY})
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${_file} OR EXISTS ${_file})
      list(APPEND _test_all_files "${_file}")
    else()
      message("The file \"${_file}\" registred by the test \"${test_name}\" does not exists")
    endif()
  endforeach()

  # add the different dependencies files (meshes, local libraries, ...)
  foreach(_dep ${_register_test_DEPENDS})
    get_target_list_of_associated_files(${_dep} _dep_ressources)
    if(_dep_ressources)
      list(APPEND _test_all_files "${_dep_ressources}")
    endif()
  endforeach()

  # add extra files to the list of files referenced by a given test
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.sh")
    list(APPEND _test_all_files "${test_name}.sh")
  endif()
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.verified")
    list(APPEND _test_all_files "${test_name}.verified")
  endif()
  if(_register_test_SCRIPT)
    list(APPEND _test_all_files "${_register_test_SCRIPT}")
  endif()

  # clean the list of all files for this test and add them in the total list
  foreach(_file ${_test_all_files})
    get_filename_component(_full ${_file} ABSOLUTE)
    file(RELATIVE_PATH __file ${PROJECT_SOURCE_DIR} ${_full})
    list(APPEND _tmp "${__file}")
  endforeach()

  foreach(_pkg ${_register_test_PACKAGE})
    package_get_name(${_pkg} _pkg_name)
    _package_add_to_variable(TESTS_FILES ${_pkg_name} ${_tmp})
  endforeach()
endfunction()
