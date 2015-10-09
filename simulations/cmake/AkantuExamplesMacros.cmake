#===============================================================================
# @file   AkantuTestAndExamples.cmake
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @date   Mon Oct 25 09:46:59 2010
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
# @section DESCRIPTION
#
#===============================================================================

#===============================================================================
macro(manage_examples et_name desc build_all label)
  string(TOUPPER ${et_name} upper_name)

  option(AKANTU_BUILD${label}${upper_name} "${desc}")
  mark_as_advanced(AKANTU_BUILD_${upper_name})

  if(${build_all})
    set(AKANTU_BUILD${label}${upper_name}_OLD
      ${AKANTU_BUILD${label}${upper_name}}
      CACHE INTERNAL "${desc}" FORCE)

    set(AKANTU_BUILD${label}${upper_name} ON
      CACHE INTERNAL "${desc}" FORCE)
  else(${build_all})
    if(DEFINED AKANTU_BUILD${label}${upper_name}_OLD)
      set(AKANTU_BUILD${label}${upper_name}
	${AKANTU_BUILD${label}${upper_name}_OLD}
	CACHE BOOL "${desc}" FORCE)

      unset(AKANTU_BUILD${label}${upper_name}_OLD
	CACHE)
    endif(DEFINED AKANTU_BUILD${label}${upper_name}_OLD)
  endif(${build_all})

  if(AKANTU_BUILD${label}${upper_name})
    add_subdirectory(${et_name})
  endif(AKANTU_BUILD${label}${upper_name})
endmacro()

#===============================================================================
# Examples
#===============================================================================
if(AKANTU_EXAMPLES)
  option(AKANTU_EXAMPLES_BUILD_ALL "Build all examples")
endif(AKANTU_EXAMPLES)

#===============================================================================
function(register_example example_name)
     set(multi_variables
    SOURCES FILES_TO_COPY DEPENDS DIRECTORIES_TO_CREATE COMPILE_OPTIONS USE
    )

  cmake_parse_arguments(_opt_pkg
    ""
    ""
    "${multi_variables}"
    ${ARGN}
    )


  set(_deps_OK TRUE)
  if(_opt_pkg_USE)
    foreach(_use ${_opt_pkg_USE})
      string(TOUPPER ${_use} _u_use)
      if(AKANTU_HAS_${_u_use})
        list(APPEND _example_INCLUDE_DIRS ${AKANTU_${_u_use}_INCLUDE_DIR})
        list(APPEND _example_LIBRARIES ${AKANTU_${_u_use}_LIBRARIES})
      else()
        message("${example_name} use ${_use} but Akantu has been compiled without this package")
        set(_deps_OK FALSE)
      endif()
    endforeach()
  endif()

  if(_deps_OK)
    include_directories(${_example_INCLUDE_DIRS})
    add_executable(${example_name} ${_opt_pkg_UNPARSED_ARGUMENTS} ${_opt_pkg_SOURCES})
    target_link_libraries(${example_name} akantu ${_example_LIBRARIES})
    if(_opt_pkg_DEPENDS)
      add_dependencies(${example_name} ${_opt_pkg_DEPENDS})
    endif()


    if(_opt_pkg_COMPILE_OPTIONS)
      set_target_properties(${test_name}
	PROPERTIES COMPILE_DEFINITIONS "${_opt_pkg_COMPILE_OPTIONS}")
    endif()

    # copy the needed files to the build folder
    if(_opt_pkg_FILES_TO_COPY)
      foreach(_file ${_opt_pkg_FILES_TO_COPY})
	file(COPY ${_file} DESTINATION .)
      endforeach()
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
macro(add_example example_name desc)
  manage_examples(${example_name} ${desc} AKANTU_EXAMPLES_BUILD_ALL _EXAMPLE_)
endmacro()
