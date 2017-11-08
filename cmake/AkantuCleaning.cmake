#===============================================================================
# @file   AkantuCleaning.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Thu Jun 1, 2017
#
# @brief  set of tools to clean the code
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
#===============================================================================

# Adding clang-format target if executable is found
find_program(CLANG_FORMAT "clang-format")
mark_as_advanced(CLANG_FORMAT)
macro(register_code_to_format)
  if(CLANG_FORMAT)
    add_custom_target(
      clang-format-all
      COMMAND ${CLANG_FORMAT}
      -i
      -style=file
      ${ARGN}
      )
  endif()
endmacro()

# Adding clang-tidy target if executable is found
find_program(CLANG_TIDY "clang-tidy")
mark_as_advanced(CLANG_TIDY)
macro(register_target_to_tidy target)
  if(CLANG_TIDY)
    option(AKANTU_CLANG_TIDY_AUTOFIX OFF)
    mark_as_advanced(AKANTU_CLANG_TIDY_AUTOFIX)

    set(_autofix_option)
    if(AKANTU_CLANG_TIDY_AUTOFIX)
      set(_autofix_option -fix)
    endif()
    get_target_property(_sources ${target} SOURCES)

    set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL
      "Enable/Disable output of compile commands during generation" FORCE)

    file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/clang-tidy)

    set(_depends)

    foreach(_src ${_sources})
      get_filename_component(_src_dir ${_src} DIRECTORY)
      file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/clang-tidy/${_src_dir})

      add_custom_command(
        OUTPUT ${PROJECT_BINARY_DIR}/clang-tidy/${_src}.yaml
        COMMAND ${CLANG_TIDY}
         -p=${PROJECT_BINARY_DIR}
         -export-fixes=${PROJECT_BINARY_DIR}/clang-tidy/${_src}.yaml
         ${_autofix_option}
        ${_src}
        COMMENT "Tidying ${_src}"
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )

      list(APPEND _depends ${PROJECT_BINARY_DIR}/clang-tidy/${_src}.yaml)
    endforeach()
    add_custom_target(clang-tidy DEPENDS ${_depends})
  endif()
endmacro()
