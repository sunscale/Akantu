#===============================================================================
# @file   AkantuTestAndExamples.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Oct 25 2010
# @date last modification: Tue Jun 24 2014
#
# @brief  macros for examples
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

#===============================================================================
function(register_example example_name)
  set(multi_variables
    SOURCES
    FILES_TO_COPY
    DEPENDS
    DIRECTORIES_TO_CREATE
    COMPILE_OPTIONS
    USE
    )

  cmake_parse_arguments(_opt_pkg
    ""
    ""
    "${multi_variables}"
    ${ARGN}
    )

  set(_deps_OK TRUE)
  if(_opt_pkg_USE)
    foreach(_pkg ${_opt_pkg_USE})
      package_is_activated(${_pkg} _activated)
      if(_activated)
        package_get_include_dir(${_pkg} _inc_dir)
        list(APPEND _example_INCLUDE_DIRS ${_inc_dir})

        package_get_libraries(${_pkg} _libraries)
        list(APPEND _example_LIBRARIES ${_libraries})

        package_get_compile_flags(${_pkg} _flags)
        list(APPEND _example_COMPILE_FLAGS "${_flags}")
      else()
        message("${example_name} use ${_pkg} but Akantu "
          " has been compiled without this package")
        set(_deps_OK FALSE)
      endif()
    endforeach()
  endif()

  if(_deps_OK)
    add_executable(${example_name}
      ${_opt_pkg_UNPARSED_ARGUMENTS} ${_opt_pkg_SOURCES})
    target_link_libraries(${example_name}
      akantu ${_example_LIBRARIES})
    target_include_directories(${example_name}
      PRIVATE ${_example_INCLUDE_DIRS})

    if(_opt_pkg_DEPENDS)
      foreach(_deps ${_opt_pkg_DEPENDS})
        get_target_property(_type ${_deps} TYPE)
        if(_type STREQUAL "SHARED_LIBRARY"
            OR _type STREQUAL "STATIC_LIBRARY")
          target_link_libraries(${example_name} ${_deps})
        else()
          add_dependencies(${example_name} ${_deps})
        endif()

      endforeach()

    endif()

    if(_opt_pkg_COMPILE_OPTIONS OR _example_COMPILE_FLAGS)
      set_target_properties(${test_name}
	PROPERTIES COMPILE_FLAGS "${_example_COMPILE_FLAGS} ${_opt_pkg_COMPILE_OPTIONS}")
    endif()

    # copy the needed files to the build folder
    if(_opt_pkg_FILES_TO_COPY)
      file(COPY ${_opt_pkg_FILES_TO_COPY} DESTINATION .)
    endif()

    # create the needed folders in the build folder
    if(_opt_pkg_DIRECTORIES_TO_CREATE)
      foreach(_dir ${_opt_pkg_DIRECTORIES_TO_CREATE})
	if(IS_ABSOLUTE ${dir})
	  file(MAKE_DIRECTORY ${_dir})
	else()
  	  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${_dir})
	endif()
      endforeach()
    endif()
  endif()
endfunction()

#===============================================================================
function(add_example et_name desc)
  string(TOUPPER ${et_name} upper_name)

  if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${et_name} AND _activated)
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
  endif()

  if(NOT (_act OR AKANTU_BUILD_ALL_EXAMPLES))
    file(RELATIVE_PATH _dir ${PROJECT_SOURCE_DIR}  ${CMAKE_CURRENT_SOURCE_DIR}/${et_name})
    list(APPEND AKANTU_TESTS_EXCLUDE_FILES /${_dir})
    set(AKANTU_TESTS_EXCLUDE_FILES ${AKANTU_TESTS_EXCLUDE_FILES} CACHE INTERNAL "")
    return()
  endif()

  option(AKANTU_BUILD_EXAMPLE_${upper_name} "${desc}")
  mark_as_advanced(AKANTU_BUILD_${upper_name})

  if(AKANTU_BUILD_ALL_EXAMPLES OR _act)
    set(AKANTU_BUILD_EXAMPLE_${upper_name}_OLD
      ${AKANTU_BUILD_EXAMPLE_${upper_name}}
      CACHE INTERNAL "${desc}" FORCE)

    set(AKANTU_BUILD_EXAMPLE_${upper_name} ${_activated}
      CACHE INTERNAL "${desc}" FORCE)
  else()
    if(DEFINED AKANTU_BUILD_EXAMPLE_${upper_name}_OLD)
      set(AKANTU_BUILD_EXAMPLE_${upper_name}
	${AKANTU_BUILD_EXAMPLE_${upper_name}_OLD}
	CACHE BOOL "${desc}" FORCE)

      unset(AKANTU_BUILD_EXAMPLE_${upper_name}_OLD
	CACHE)
    endif(DEFINED AKANTU_BUILD_EXAMPLE_${upper_name}_OLD)
  endif()

  if(AKANTU_BUILD_EXAMPLE_${upper_name})
    add_subdirectory(${et_name})
  endif(AKANTU_BUILD_EXAMPLE_${upper_name})
endfunction()
