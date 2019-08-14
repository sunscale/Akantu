#===============================================================================
# @file   AkantuMacros.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Oct 22 2010
# @date last modification: Tue Jan 19 2016
#
# @brief  Set of macros used by akantu cmake files
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

#===============================================================================
function(set_third_party_shared_libirary_name _var _lib)
  set(${_var}
    ${PROJECT_BINARY_DIR}/third-party/lib/${CMAKE_SHARED_LIBRARY_PREFIX}${_lib}${CMAKE_SHARED_LIBRARY_SUFFIX}
    CACHE FILEPATH "" FORCE)
endfunction()

# ==============================================================================
function(_add_file_to_copy target file)
  get_filename_component(_file_name_we ${file} NAME_WE)
  get_filename_component(_file_name ${file} NAME)
  get_filename_component(_file_path ${file}
    ABSOLUTE BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

  set(copy_target copy_${_file_name_we}_${target})
  add_custom_target(${copy_target}
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${_file_name})
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${_file_name}
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${file}
		${CMAKE_CURRENT_BINARY_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS ${_file_path}
    COMMENT "Copying file ${_file_name} for the target ${target}"
    )
  add_dependencies(${target} ${copy_target})
endfunction()

# ==============================================================================
function(get_target_list_of_associated_files tgt files)
  if(TARGET ${tgt})
    get_target_property(_type ${tgt} TYPE)
  else()
    set(_type ${tgt}-NOTFOUND)
  endif()

  if(_type STREQUAL "SHARED_LIBRARY"
      OR _type STREQUAL "STATIC_LIBRARY"
      OR _type STREQUAL "MODULE_LIBRARY"
      OR _type STREQUAL "EXECUTABLE")
    get_target_property(_srcs ${tgt} SOURCES)
    set(_dep_ressources)
    foreach(_file ${_srcs})
      list(APPEND _dep_ressources ${CMAKE_CURRENT_SOURCE_DIR}/${_file})
    endforeach()
  elseif(_type)
    get_target_property(_dep_ressources ${tgt} RESSOURCES)
  endif()

  set(${files} ${_dep_ressources} PARENT_SCOPE)
endfunction()

#===============================================================================
# Generate the list of currently loaded materials
function(generate_material_list)
  message(STATUS "Determining the list of recognized materials...")

  package_get_all_include_directories(
    AKANTU_LIBRARY_INCLUDE_DIRS
    )

  package_get_all_external_informations(
    PRIVATE_INCLUDE AKANTU_PRIVATE_EXTERNAL_INCLUDE_DIR
    INTERFACE_INCLUDE AKANTU_INTERFACE_EXTERNAL_INCLUDE_DIR
    LIBRARIES AKANTU_EXTERNAL_LIBRARIES
    )

  set(_include_dirs
    ${AKANTU_INCLUDE_DIRS}
    ${AKANTU_PRIVATE_EXTERNAL_INCLUDE_DIR}
    ${AKANTU_INTERFACE_EXTERNAL_INCLUDE_DIR}
    )

  try_run(_material_list_run _material_list_compile
    ${CMAKE_BINARY_DIR}
    ${PROJECT_SOURCE_DIR}/cmake/material_lister.cc
    CMAKE_FLAGS "-DINCLUDE_DIRECTORIES:STRING=${_include_dirs}" "-DCMAKE_CXX_STANDARD=14"
    COMPILE_DEFINITIONS "-DAKANTU_CMAKE_LIST_MATERIALS"
    COMPILE_OUTPUT_VARIABLE _compile_results
    RUN_OUTPUT_VARIABLE _result_material_list)

  if(_material_list_compile AND "${_material_list_run}" EQUAL 0)
    message(STATUS "Materials included in Akantu:")
    string(REPLACE "\n" ";" _material_list "${_result_material_list}")
    foreach(_mat ${_material_list})
      string(REPLACE ":" ";" _mat_key "${_mat}")
      list(GET _mat_key 0 _key)
      list(GET _mat_key 1 _class)
      list(LENGTH _mat_key _l)

      if("${_l}" GREATER 2)
        list(REMOVE_AT _mat_key 0 1)
        set(_opt " -- options: [")
        foreach(_o ${_mat_key})
          set(_opt "${_opt} ${_o}")
        endforeach()
        set(_opt "${_opt} ]")
      else()
        set(_opt "")
      endif()

      message(STATUS "   ${_class} -- key: ${_key}${_opt}")
    endforeach()
  else()
    message(STATUS "Could not determine the list of materials.")
    message("${_compile_results}")
  endif()
endfunction()

#===============================================================================
# Declare the options for the types and defines the approriate typedefs
function(declare_akantu_types)
  set(AKANTU_TYPE_FLOAT "double (64bit)" CACHE STRING "Precision force floating point types")
  mark_as_advanced(AKANTU_TYPE_FLOAT)
  set_property(CACHE AKANTU_TYPE_FLOAT PROPERTY STRINGS
    "quadruple (128bit)"
    "double (64bit)"
    "float (32bit)"
    )

  set(AKANTU_TYPE_INTEGER "int (32bit)" CACHE STRING "Size of the integer types")
  mark_as_advanced(AKANTU_TYPE_INTEGER)
  set_property(CACHE AKANTU_TYPE_INTEGER PROPERTY STRINGS
    "int (32bit)"
    "long int (64bit)"
    )

  include(CheckTypeSize)

  # ----------------------------------------------------------------------------
  # Floating point types
  # ----------------------------------------------------------------------------
  if(AKANTU_TYPE_FLOAT STREQUAL "float (32bit)")
    set(AKANTU_FLOAT_TYPE "float" CACHE INTERNAL "")
    set(AKANTU_FLOAT_SIZE 4 CACHE INTERNAL "")
  elseif(AKANTU_TYPE_FLOAT STREQUAL "double (64bit)")
    set(AKANTU_FLOAT_TYPE "double" CACHE INTERNAL "")
    set(AKANTU_FLOAT_SIZE 8 CACHE INTERNAL "")
  elseif(AKANTU_TYPE_FLOAT STREQUAL "quadruple (128bit)")
    check_type_size("long double" LONG_DOUBLE)
    if(HAVE_LONG_DOUBLE)
      set(AKANTU_FLOAT_TYPE "long double" CACHE INTERNAL "")
      set(AKANTU_FLOAT_SIZE 16 CACHE INTERNAL "")
      message("This feature is not tested and will most probably not compile")
    else()
      message(FATAL_ERROR "The type long double is not defined on your system")
    endif()
  else()
    message(FATAL_ERROR "The float type is not defined")
  endif()

  include(CheckIncludeFileCXX)
  include(CheckCXXSourceCompiles)

  # ----------------------------------------------------------------------------
  # Integer types
  # ----------------------------------------------------------------------------
  check_include_file_cxx(cstdint HAVE_CSTDINT)
  if(NOT HAVE_CSTDINT)
    check_include_file_cxx(stdint.h HAVE_STDINT_H)
    if(HAVE_STDINT_H)
      list(APPEND _int_include stdint.h)
    endif()
  else()
    list(APPEND _int_include cstdint)
  endif()


  check_include_file_cxx(cstddef HAVE_CSTDDEF)
  if(NOT HAVE_CSTDINT)
    check_include_file_cxx(stddef.h HAVE_STDDEF_H)
    if(HAVE_STDINT_H)
      list(APPEND _int_include stddef.h)
    endif()
  else()
    list(APPEND _int_include cstddef)
  endif()

  if(AKANTU_TYPE_INTEGER STREQUAL "int (32bit)")
    set(AKANTU_INTEGER_SIZE 4 CACHE INTERNAL "")
    check_type_size("int" INT)
    if(INT EQUAL 4)
      set(AKANTU_SIGNED_INTEGER_TYPE "int" CACHE INTERNAL "")
      set(AKANTU_UNSIGNED_INTEGER_TYPE "unsigned int" CACHE INTERNAL "")
    else()
      check_type_size("int32_t" INT32_T LANGUAGE CXX)
      if(HAVE_INT32_T)
        set(AKANTU_SIGNED_INTEGER_TYPE "int32_t" CACHE INTERNAL "")
        set(AKANTU_UNSIGNED_INTEGER_TYPE "uint32_t" CACHE INTERNAL "")
        list(APPEND _extra_includes ${_int_include})
      endif()
    endif()
  elseif(AKANTU_TYPE_INTEGER STREQUAL "long int (64bit)")
    set(AKANTU_INTEGER_SIZE 8 CACHE INTERNAL "")
    check_type_size("long int" LONG_INT)
    if(LONG_INT EQUAL 8)
      set(AKANTU_SIGNED_INTEGER_TYPE "long int" CACHE INTERNAL "")
      set(AKANTU_UNSIGNED_INTEGER_TYPE "unsigned long int" CACHE INTERNAL "")
    else()
      check_type_size("long long int" LONG_LONG_INT)
      if(HAVE_LONG_LONG_INT AND LONG_LONG_INT EQUAL 8)
        set(AKANTU_SIGNED_INTEGER_TYPE "long long int" CACHE INTERNAL "")
        set(AKANTU_UNSIGNED_INTEGER_TYPE "unsigned long long int" CACHE INTERNAL "")
      else()
        check_type_size("int64_t" INT64_T)
        if(HAVE_INT64_T)
          set(AKANTU_SIGNED_INTEGER_TYPE "int64_t" CACHE INTERNAL "")
          set(AKANTU_UNSIGNED_INTEGER_TYPE "uint64_t" CACHE INTERNAL "")
          list(APPEND _extra_includes ${_int_include})
        endif()
      endif()
    endif()
  else()
    message(FATAL_ERROR "The integer type is not defined")
  endif()

  # ----------------------------------------------------------------------------
  # includes
  # ----------------------------------------------------------------------------
  foreach(_inc ${_extra_includes})
    set(_incs "#include <${_inc}>\n${_incs}")
  endforeach()
  set(AKANTU_TYPES_EXTRA_INCLUDES ${_incs} CACHE INTERNAL "")
endfunction()


function(mask_package_options prefix)
  get_property(_list DIRECTORY PROPERTY VARIABLES)
  foreach(_var ${_list})
    if("${_var}" MATCHES "^${prefix}.*")
      mark_as_advanced(${_var})
    endif()
  endforeach()
endfunction()
