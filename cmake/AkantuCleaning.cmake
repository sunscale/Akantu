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
find_program(CLANG_FORMAT_EXECUTABLE "clang-format")
# Adding clang-tidy target if executable is found
find_program(CLANG_TIDY_EXECUTABLE "clang-tidy")
find_program(RUN_CLANG_TIDY_EXECUTABLE "run-clang-tidy")

mark_as_advanced(
  CLANG_FORMAT_EXECUTABLE
  CLANG_TIDY_EXECUTABLE
  RUN_CLANG_TIDY_EXECUTABLE
  )

function(register_code_to_format)
  if(NOT CLANG_FORMAT_EXECUTABLE)
  endif()

  add_custom_target(
    clang-format-all
    COMMAND ${CLANG_FORMAT_EXECUTABLE}
    -style=file
    --output-replacements-xml
    ${ARGN} || /bin/true
    )
endfunction()

function(register_tidy_all directory)
  if(NOT RUN_CLANG_TIDY_EXECUTABLE)
    return()
  endif()

  add_custom_target(
    clang-tidy-all
    COMMAND ${RUN_CLANG_TIDY_EXECUTABLE}
    ${ARGN}
    )
endfunction()

function(register_target_to_tidy target)
  if(NOT CLANG_TIDY_EXECUTABLE)
    return()
  endif()

  option(AKANTU_CLANG_TIDY_AUTOFIX OFF)
  mark_as_advanced(AKANTU_CLANG_TIDY_AUTOFIX)

  set(_autofix_option)
  if(AKANTU_CLANG_TIDY_AUTOFIX)
    set(_autofix_option -fix)
  endif()
  get_target_property(_sources ${target} SOURCES)

  file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/clang-tidy)

  set(_depends)

  foreach(_src ${_sources})
    get_filename_component(_src_dir ${_src} DIRECTORY)
    file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/clang-tidy/${_src_dir})

    add_custom_command(
      OUTPUT ${PROJECT_BINARY_DIR}/clang-tidy/${_src}.yaml
      COMMAND ${CLANG_TIDY_EXECUTABLE}
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
endfunction()
