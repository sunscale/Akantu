#===============================================================================
# @file   AkantuExampleMacros.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Jan 18 2016
# @date last modification: Fri Jan 22 2016
#
# @brief  macros for examples
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
# @section DESCRIPTION
#
#===============================================================================

include(AkantuSimulationMacros)

# ==============================================================================
function(register_example example_name)
  _add_akantu_simulation(${example_name} ${ARGN} LIST_FILES _example_files)
  if(DEFINED _add_examples_pkg)
    package_add_to_variable(EXAMPLES_FILES ${_add_examples_pkg} ${_example_files})
  endif()

  if(AKANTU_TEST_EXAMPLES)
    cmake_parse_arguments(_example
      "PYTHON;PARALLEL"
      ""
      "SCRIPT"
      ${ARGN}
      )

    if(_example_PARALLEL)
      set(_exe ${MPIEXEC})
      if(NOT _exe)
	set(_exe ${MPIEXEC_EXECUTABLE})
      endif()
      set(_parallel_runner ${_exe} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 2)
    endif()
    
    if(NOT _example_SCRIPT)
      add_test(NAME ${example_name}-test
	COMMAND ${_parallel_runner} $<TARGET_FILE:${example_name}>
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    elseif(_example_SCRIPT)
      if(_example_PYTHON)
	add_test(NAME ${example_name}-test
	  COMMAND ${_parallel_runner} ${PYTHON_EXECUTABLE} ${_example_SCRIPT}
	  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
      else()
	set(_python_path ENV{PYTHON_PATH})
	if (NOT _python_path MATCHES ${PROJECT_BINARY_DIR}/python)
	  set(ENV{PYTHON_PATH} "${_python_path}:${PROJECT_BINARY_DIR}/python")
	endif()
	add_test(NAME ${example_name}-test
	  COMMAND ${_parallel_runner} ${_example_SCRIPT}
	  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	  )
      endif()
    endif()
  endif()
endfunction()

# ==============================================================================
function(add_example et_name desc)
  string(TOUPPER ${et_name} upper_name)

  if(NOT _build_all_ex)
    option(AKANTU_BUILD_ALL_EXAMPLES "Activate all examples" OFF)
    set( _build_all_ex TRUE)
  endif()

  option(AKANTU_BUILD_EXAMPLES_${upper_name} ${desc} OFF)

  if(AKANTU_BUILD_ALL_EXAMPLES)
    mark_as_advanced(FORCE AKANTU_BUILD_EXAMPLES_${upper_name})
  else()
    mark_as_advanced(CLEAR AKANTU_BUILD_EXAMPLES_${upper_name})
  endif()

  if(AKANTU_BUILD_EXAMPLES_${upper_name} OR AKANTU_BUILD_ALL_EXAMPLES)

    if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${et_name})
      message(FATAL_ERROR "The folder ${CMAKE_CURRENT_SOURCE_DIR}/${et_name} "
	"that you try to register as an example sub-folder, does not exists.")
    endif()

    cmake_parse_arguments(_manage_example
      ""
      ""
      "PACKAGE"
      ${ARGN}
      )

    if(_manage_example_PACKAGE)
      set(_act TRUE)
      foreach(_pkg ${_manage_example_PACKAGE})
	package_is_activated(${_pkg} _activated)
	if(NOT _activated)
          set(_act FALSE)
	endif()
      endforeach()
    else()
      message(SEND_ERROR "Examples should be associated to a package")
    endif()

    if(_act)
      if(DEFINED _add_examples_pkg)
	set(_save_add_examples_pkg ${_add_examples_pkg})
      endif()
      list(GET _manage_example_PACKAGE 0 _pkg)
      set(_add_examples_pkg ${_pkg})

      add_subdirectory(${et_name})

      unset(_add_examples_pkg)
      if(DEFINED _save_add_examples_pkg)
	set(_add_examples_pkg ${_save_add_examples_pkg})
	unset(_save_add_examples_pkg)
      endif()
    endif()
  endif()
endfunction()
