#===============================================================================
# @file   CMakeDetermineCompiler.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Sun Oct 19 2014
#
# @brief  
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
#===============================================================================

# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

macro(_cmake_find_compiler lang)
  # Use already-enabled languages for reference.
  get_property(_languages GLOBAL PROPERTY ENABLED_LANGUAGES)
  list(REMOVE_ITEM _languages "${lang}")

  if(CMAKE_${lang}_COMPILER_INIT)
    # Search only for the specified compiler.
    set(CMAKE_${lang}_COMPILER_LIST "${CMAKE_${lang}_COMPILER_INIT}")
  else()
    # Re-order the compiler list with preferred vendors first.
    set(_${lang}_COMPILER_LIST "${CMAKE_${lang}_COMPILER_LIST}")
    set(CMAKE_${lang}_COMPILER_LIST "")
    # Prefer vendors of compilers from reference languages.
    foreach(l ${_languages})
      list(APPEND CMAKE_${lang}_COMPILER_LIST
        ${_${lang}_COMPILER_NAMES_${CMAKE_${l}_COMPILER_ID}})
    endforeach()
    # Prefer vendors based on the platform.
    list(APPEND CMAKE_${lang}_COMPILER_LIST ${CMAKE_${lang}_COMPILER_NAMES})
    # Append the rest of the list and remove duplicates.
    list(APPEND CMAKE_${lang}_COMPILER_LIST ${_${lang}_COMPILER_LIST})
    unset(_${lang}_COMPILER_LIST)
    list(REMOVE_DUPLICATES CMAKE_${lang}_COMPILER_LIST)
    if(CMAKE_${lang}_COMPILER_EXCLUDE)
      list(REMOVE_ITEM CMAKE_${lang}_COMPILER_LIST
        ${CMAKE_${lang}_COMPILER_EXCLUDE})
    endif()
  endif()

  # Look for directories containing compilers of reference languages.
  set(_${lang}_COMPILER_HINTS)
  foreach(l ${_languages})
    if(CMAKE_${l}_COMPILER AND IS_ABSOLUTE "${CMAKE_${l}_COMPILER}")
      get_filename_component(_hint "${CMAKE_${l}_COMPILER}" PATH)
      if(IS_DIRECTORY "${_hint}")
        list(APPEND _${lang}_COMPILER_HINTS "${_hint}")
      endif()
      unset(_hint)
    endif()
  endforeach()

  # Find the compiler.
  if(_${lang}_COMPILER_HINTS)
    # Prefer directories containing compilers of reference languages.
    list(REMOVE_DUPLICATES _${lang}_COMPILER_HINTS)
    find_program(CMAKE_${lang}_COMPILER
      NAMES ${CMAKE_${lang}_COMPILER_LIST}
      PATHS ${_${lang}_COMPILER_HINTS}
      NO_DEFAULT_PATH
      DOC "${lang} compiler")
  endif()
  find_program(CMAKE_${lang}_COMPILER NAMES ${CMAKE_${lang}_COMPILER_LIST} DOC "${lang} compiler")
  if(CMAKE_${lang}_COMPILER_INIT AND NOT CMAKE_${lang}_COMPILER)
    set(CMAKE_${lang}_COMPILER "${CMAKE_${lang}_COMPILER_INIT}" CACHE FILEPATH "${lang} compiler" FORCE)
  endif()
  unset(_${lang}_COMPILER_HINTS)
  unset(_languages)

  # Look for a make tool provided by Xcode
  if(CMAKE_${lang}_COMPILER STREQUAL "CMAKE_${lang}_COMPILER-NOTFOUND" AND CMAKE_HOST_APPLE)
    foreach(comp ${CMAKE_${lang}_COMPILER_LIST})
      execute_process(COMMAND xcrun --find ${comp}
        OUTPUT_VARIABLE _xcrun_out OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_VARIABLE _xcrun_err)
      if(_xcrun_out)
        set_property(CACHE CMAKE_${lang}_COMPILER PROPERTY VALUE "${_xcrun_out}")
        break()
      endif()
    endforeach()
  endif()
endmacro()
