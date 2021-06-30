#===============================================================================
# @file   CMakePackagesSystem.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Wed Nov 05 2014
# @date last modification: Wed Jan 20 2016
#
# @brief  Set of macros used by akantu to handle the package system
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

#[=======================================================================[.rst:
#CMakePackagesSystem
#-------------------
#
#This package defines multiple function to handle packages. This packages can
#be of two kinds regular ones and extra_packages (ex: in akantu the LGPL part
#is regular packages and extra packages are on Propetary license)
#
#Package are loaded with the help of the command:
#
#.. command:: package_list_packages
#
#     package_list_packages(<regular_package_folder>
#       [ EXTRA_PACKAGE_FOLDER <extra_package_folder> ]
#       [ SOURCE_FOLDER <source_folder>]
#       [ TEST_FOLDER <test_folder> ]
#       [ MANUAL_FOLDER <manual_folder> ]
#       )
#
#     This command will look for packages name like ``<regular_package_folder>/<package>.cmake``
#     OR ``<extra_package_folder>/<package>/package.cmake``
#
#A package is a cmake script that should contain at list the declaration of a
#package
#
#.. command:: package_declare
#
#     package_declare(<package real name>
#       [EXTERNAL] [META] [ADVANCED] [NOT_OPTIONAL]
#       [DESCRIPTION <description>] [DEFAULT <default_value>]
#       [DEPENDS <pkg> ...]
#       [BOOST_COMPONENTS <pkg> ...]
#       [EXTRA_PACKAGE_OPTIONS <opt> ...]
#       [COMPILE_FLAGS <lang> <flags>]
#       [SYSTEM <ON|OFF|AUTO> [ <script_to_compile> ]]
#       [FEATURES_PUBLIC <feature> ...]
#       [FEATURES_PRIVATE <feature> ...]
#       [EXCLUDE_FROM_ALL]
#       )
#
#.. command:: package_declare_sources
#
#     It can also declare multiple informations:
#     source files:
#
#     package_declare_sources(<package real name>
#       <src1> <src2> ... <srcn>)
#
#.. command:: package_declare_documentation
#
#     a LaTeX documentation
#     package_declare_documentation(<package real name>
#       <line1> <line2> ...<linen>)
#
#.. command:: package_declare_documentation_files
#
#     LaTeX documentation files
#     package_declare_documentation_files(<package real name>
#       <file1> <file2> ... <filen>)
#
#Different function can also be retrieved from the package system by using the
#different accessors
#
#.. command:: package_get_name
#     package_get_name(<pkg> <retval>)
#
#.. command:: package_get_real_name
#    package_get_real_name(<pkg> <retval>)
#
#.. command:: package_get_option_name
#     package_get_option_name(<pkg> <retval>)
#
#.. command:: package_use_system
#     package_use_system(<pkg> <retval>)
#
#.. command:: package_get_nature
#     package_get_nature(<pkg> <retval>)
#
#.. command:: package_get_description
#     package_get_description(<pkg> <retval>)
#
#.. command:: package_get_filename
#     package_get_filename(<pkg> <retval>)
#
#.. command:: package_get_sources_folder
#     package_get_sources_folder(<pkg> <retval>)
#.. command:: package_get_tests_folder
#     package_get_tests_folder(<pkg> <retval>)
#.. command:: package_get_manual_folder
#     package_get_manual_folder(<pkg> <retval>)
#
#.. command:: package_get_find_package_extra_options
#     package_get_find_package_extra_options(<pkg> <retval>)
#
#.. command:: package_get_compile_flags
#     package_get_compile_flags(<pkg> <lang> <retval>)
#.. command:: package_set_compile_flags
#     package_set_compile_flags(<pkg> <lang> <flag1> <flag2> ... <flagn>)
#
#.. command:: package_get_include_dir
#     package_get_include_dir(<pkg> <retval>)
#.. command:: package_set_include_dir
#     package_set_include_dir(<pkg> <inc1> <inc2> ... <incn>)
#.. command:: package_add_include_dir
#     package_add_include_dir(<pkg> <inc1> <inc2> ... <incn>)
#
#.. command:: package_get_libraries
#     package_get_libraries(<pkg> <retval>)
#.. command:: package_set_libraries
#     package_set_libraries(<pkg> <lib1> <lib2> ... <libn>)
#
#.. command:: package_add_extra_dependency
#     package_add_extra_dependency(pkg <dep1> <dep2> ... <depn>)
#.. command:: package_rm_extra_dependency
#     package_rm_extra_dependency(<pkg> <dep>)
#.. command:: package_get_extra_dependencies
#     package_get_extra_dependencies(<pkg> <retval>)
#
#.. command:: package_is_activated
#     package_is_activated(<pkg> <retval>)
#.. command:: package_is_deactivated
#     package_is_deactivated(<pkg> <retval>)
#
#.. command:: package_get_dependencies
#     package_get_dependencies(<pkg> <PRIVATE|INTERFACE> <retval>)
#.. command:: package_add_dependencies
#     package_add_dependencies(<pkg> <PRIVATE|INTERFACE> <dep1> <dep2> ... <depn>)
#     package_remove_dependencies(<pkg> <dep1> <dep2> ... <depn>)
#     package_remove_dependency(<pkg> <dep>)
#
#.. command:: package_on_enabled_script
#     package_on_enabled_script(<pkg> <script>)
#
#.. command:: package_get_all_source_files
#     package_get_all_source_files(<srcs> <public_headers> <private_headers>)
#.. command:: package_get_all_include_directories
#     package_get_all_include_directories(<inc_dirs>)
#.. command:: package_get_all_external_informations
#     package_get_all_external_informations(<include_dir> <libraries>)
#.. command:: package_get_all_definitions
#     package_get_all_definitions(<definitions>)
#.. command:: package_get_all_extra_dependencies
#     package_get_all_extra_dependencies(<dependencies>)
#.. command:: package_get_all_test_folders
#     package_get_all_test_folders(<test_dirs>)
#.. command:: package_get_all_documentation_files
#     package_get_all_documentation_files(<doc_files>)
#.. command:: package_get_all_activated_packages
#     package_get_all_activated_packages(<activated_list>)
#.. command:: package_get_all_deactivated_packages
#     package_get_all_deactivated_packages(<deactivated_list>)
#.. command:: package_get_all_packages
#     package_get_all_packages(<packages_list>)
#.. command:: package_get_all_features_public
#     package_get_all_features_public(<features>)
#.. command:: package_get_all_features_private
#     package_get_all_features_private(<features>)
#
#
#     .. command:: package_set_package_system_dependency
#
#     package_set_package_system_dependency(<pkg> <system> <dep1>
#                                           <dep2> ... <depn>)
#
#                                       .. command:: package_get_package_system_dependency
#
#     package_get_package_system_dependency(<pkg> <var>)
#
#
#]=======================================================================]
if (DEFINED CMAKE_PACKAGES_SYSTEM_LOADED)
  return()
endif()
set(CMAKE_PACKAGES_SYSTEM_LOADED TRUE)


include(CMakeParseArguments)

#===============================================================================
# Package Management
#===============================================================================
if(__CMAKE_PACKAGES_SYSTEM)
  return()
endif()
set(__CMAKE_PACKAGES_SYSTEM TRUE)

if(CMAKE_VERSION VERSION_GREATER 3.1.2)
  cmake_policy(SET CMP0054 NEW)
endif()

#===============================================================================
option(AUTO_MOVE_UNKNOWN_FILES
  "Give to cmake the permission to move the unregistered files to the ${PROJECT_SOURCE_DIR}/tmp directory" FALSE)
mark_as_advanced(AUTO_MOVE_UNKNOWN_FILES)


include(CMakePackagesSystemGlobalFunctions)
include(CMakePackagesSystemPrivateFunctions)

# ==============================================================================
# "Public" Accessors
# ==============================================================================
# ------------------------------------------------------------------------------
# Package name
# ------------------------------------------------------------------------------
function(package_get_name pkg pkg_name)
  string(TOUPPER ${PROJECT_NAME} _project)
  string(REPLACE "-" "_" _str_pkg "${pkg}")
  string(TOUPPER ${_str_pkg} _u_package)
  set(${pkg_name} ${_project}_PKG_${_u_package} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Real name
# ------------------------------------------------------------------------------
function(package_get_real_name pkg ret)
  package_get_name(${pkg} _pkg_name)
  _package_get_real_name(${_pkg_name} _tmp)
  set(${ret} ${_tmp} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Option name
# ------------------------------------------------------------------------------
function(package_get_option_name pkg ret)
  package_get_name(${pkg} _pkg_name)
  _package_get_option_name(${_pkg_name} _tmp)
  set(${ret} ${_tmp} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Set if system package or compile external lib
# ------------------------------------------------------------------------------
function(package_use_system pkg ret)
  package_get_name(${pkg} _pkg_name)
  _package_use_system(${_pkg_name} _tmp)
  set(${ret} ${_tmp} PARENT_SCOPE)
endfunction()

function(package_add_third_party_script_variable pkg var)
  package_get_name(${pkg} _pkg_name)
  _package_add_third_party_script_variable(${_pkg_name} ${var} ${ARGN})
  set(${var} ${ARGN} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Add package's targets to the export list
# ------------------------------------------------------------------------------
function(package_add_to_export_list pkg)
  package_get_name(${pkg} _pkg_name)
  _package_add_to_export_list(${_pkg_name} ${ARGN})
endfunction()

# ------------------------------------------------------------------------------
# Removes packages's targets from export list
# ------------------------------------------------------------------------------
function(package_remove_from_export_list pkg)
  package_get_name(${pkg} _pkg_name)
  _package_remove_from_export_list(${_pkg_name} ${ARGN})
endfunction()

# ------------------------------------------------------------------------------
# Nature
# ------------------------------------------------------------------------------
function(package_get_nature pkg ret)
  package_get_name(${pkg} _pkg_name)
  _package_get_nature(${_pkg_name} _tmp)
  set(${ret} ${_tmp} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Description
# ------------------------------------------------------------------------------
function(package_get_description pkg ret)
  package_get_name(${pkg} _pkg_name)
  _package_get_description(${_pkg_name} _tmp)
  set(${ret} ${_tmp} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Package file name
# ------------------------------------------------------------------------------
function(package_get_filename pkg ret)
  package_get_name(${pkg} _pkg_name)
  _package_get_filename(${_pkg_name} _tmp)
  set(${ret} ${_tmp} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Source files
# ------------------------------------------------------------------------------
function(package_get_source_files pkg ret_srcs ret_pub ret_priv)
  package_get_name(${pkg} _pkg_name)
  _package_get_source_files(${_pkg_name} _tmp_srcs _tmp_pub _tmp_priv)
  set(${ret_srcs} ${_tmp_srcs} PARENT_SCOPE)
  set(${ret_pub}  ${_tmp_pub}  PARENT_SCOPE)
  set(${ret_priv} ${_tmp_pric} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Source folder
# ------------------------------------------------------------------------------
function(package_get_sources_folder pkg ret)
  package_get_name(${pkg} _pkg_name)
  _package_get_sources_folder(${_pkg_name} _tmp)
  set(${ret} ${_tmp} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Test folder
# ------------------------------------------------------------------------------
function(package_get_tests_folder pkg ret)
  package_get_name(${pkg} _pkg_name)
  _package_get_tests_folder(${_pkg_name} _tmp)
  set(${ret} ${_tmp} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Manual folder
# ------------------------------------------------------------------------------
function(package_get_manual_folder pkg ret)
  package_get_name(${pkg} _pkg_name)
  _package_get_manual_folder(${_pkg_name} _tmp)
  set(${ret} ${_tmp} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Extra option for the find_package
# ------------------------------------------------------------------------------
function(package_get_find_package_extra_options pkg ret)
  package_get_name(${pkg} _pkg_name)
  _package_get_find_package_extra_options(${_pkg_name} _tmp)
  set(${ret} ${_tmp} PARENT_SCOPE)
endfunction()

function(package_set_find_package_extra_options pkg)
  package_get_name(${pkg} _pkg_name)
  _package_set_find_package_extra_options(${_pkg_name} ${ARGN})
endfunction()

# ------------------------------------------------------------------------------
# Compilation flags
# ------------------------------------------------------------------------------
function(package_get_compile_flags pkg lang ret)
  package_get_name(${pkg} _pkg_name)
  _package_get_compile_flags(${_pkg_name} ${lang} _tmp)
  set(${ret} "${_tmp}" PARENT_SCOPE)
endfunction()

function(package_set_compile_flags pkg lang)
  package_get_name(${pkg} _pkg_name)
  _package_set_compile_flags(${_pkg_name} ${lang} ${ARGN})
endfunction()

function(package_unset_compile_flags pkg lang)
  package_get_name(${pkg} _pkg_name)
  _package_unset_compile_flags(${_pkg_name} ${lang})
endfunction()

# ------------------------------------------------------------------------------
# Include dir
# ------------------------------------------------------------------------------
function(package_get_include_dir pkg ret)
  package_get_name(${pkg} _pkg_name)
  _package_get_include_dir(${_pkg_name} _tmp)
  set(${ret} ${_tmp} PARENT_SCOPE)
endfunction()

function(package_set_include_dir pkg)
  package_get_name(${pkg} _pkg_name)
  _package_set_include_dir(${_pkg_name} ${ARGN})
endfunction()

function(package_add_include_dir pkg)
  package_get_name(${pkg} _pkg_name)
  _package_add_include_dir(${_pkg_name} ${ARGN})
endfunction()

# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------
function(package_get_libraries pkg ret)
  package_get_name(${pkg} _pkg_name)
  _package_get_libraries(${_pkg_name} _tmp)
  set(${ret} ${_tmp} PARENT_SCOPE)
endfunction()

function(package_set_libraries pkg)
  package_get_name(${pkg} _pkg_name)
  _package_set_libraries(${_pkg_name} ${ARGN})
endfunction()

# ------------------------------------------------------------------------------
# Extra dependencies like custom commands of ExternalProject
# ------------------------------------------------------------------------------
function(package_add_extra_dependency pkg)
  package_get_name(${pkg} _pkg_name)
  _package_add_extra_dependency(${_pkg_name} ${ARGN})
endfunction()

function(package_rm_extra_dependency pkg dep)
  package_get_name(${pkg} _pkg_name)
  _package_rm_extra_dependency(${_pkg_name} ${dep})
endfunction()

function(package_get_extra_dependencies pkg ret)
  package_get_name(${pkg} _pkg_name)
  _package_get_extra_dependencies(${_pkg_name} _tmp)
  set(${ret} ${_tmp} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Activate/deactivate
# ------------------------------------------------------------------------------
function(package_is_activated pkg ret)
  package_get_name(${pkg} _pkg_name)
  _package_is_activated(${_pkg_name} _tmp)
  set(${ret} ${_tmp} PARENT_SCOPE)
endfunction()

function(package_is_deactivated pkg ret)
  package_get_name(${pkg} _pkg_name)
  _package_is_deactivated(${_pkg_name} _tmp)
  set(${ret} ${_tmp} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Direct dependencies
# ------------------------------------------------------------------------------
function(package_get_dependencies pkg type ret)
  package_get_name(${pkg} _pkg_name)
  _package_get_dependencies(${_pkg_name} ${type} _tmp_name)
  _package_get_real_name(${_tmp_name} _tmp)
  set(${ret} ${_tmp} PARENT_SCOPE)
endfunction()

function(package_add_dependencies pkg type)
  package_get_name(${pkg} _pkg_name)
  foreach(_dep ${ARGN})
    package_get_name(${_dep} _dep_pkg_name)
    list(APPEND _tmp_deps ${_dep_pkg_name})
  endforeach()

  _package_add_dependencies(${_pkg_name} ${type} ${_tmp_deps})
endfunction()

function(package_remove_dependencies pkg type)
  foreach(_dep ${ARGN})
    package_remove_dependency(${pkg} _dep)
  endforeach()
endfunction()

function(package_remove_dependency pkg dep)
  package_get_name(${pkg} _pkg_name)
  package_get_name(${dep} _dep_pkg_name)
  _package_remove_dependency(${_pkg_name} PRIVATE ${_dep_pkg_name})
  _package_remove_dependency(${_pkg_name} INTERFACE ${_dep_pkg_name})
endfunction()

# ------------------------------------------------------------------------------
# Documentation related functions
# ------------------------------------------------------------------------------
function(package_declare_documentation pkg)
  package_get_name(${pkg} _pkg_name)
  _package_set_documentation(${_pkg_name} ${ARGN})
endfunction()

function(package_declare_documentation_files pkg)
  package_get_name(${pkg} _pkg_name)
  _package_set_documentation_files(${_pkg_name} ${ARGN})
endfunction()

# ------------------------------------------------------------------------------
# Set any user variables needed
# ------------------------------------------------------------------------------
function(package_set_variable variable pkg)
  package_get_name(${pkg} _pkg_name)
  _package_set_variable(${variable} ${_pkg_name} ${ARGN})
endfunction()

function(package_add_to_variable variable pkg)
  package_get_name(${pkg} _pkg_name)
  _package_add_to_variable(${variable} ${_pkg_name} ${ARGN})
endfunction()

function(package_get_variable variable pkg value)
  package_get_name(${pkg} _pkg_name)
  _package_get_variable(${variable} ${_pkg_name} _value_tmp)
  if(_value_tmp)
    set(${value} ${_value_tmp} PARENT_SCOPE)
  else()
    unset(${value} PARENT_SCOPE)
  endif()
endfunction()

# ------------------------------------------------------------------------------
# Exteral package system as apt rpm dependencies
# ------------------------------------------------------------------------------
function(package_set_package_system_dependency pkg system)
  package_get_name(${pkg} _pkg_name)
  _package_set_package_system_dependency(${_pkg_name} ${system} ${ARGN})
endfunction()

function(package_get_package_system_dependency pkg system var)
  package_get_name(${pkg} _pkg_name)
  _package_set_package_system_dependency(${_pkg_name} ${sytem} _tmp)
  set(${var} ${_tmp} PARENT_SCOPE)
endfunction()
# ------------------------------------------------------------------------------

# ==============================================================================
# Global accessors
# ==============================================================================
# ------------------------------------------------------------------------------
# get the list of source files
# ------------------------------------------------------------------------------
function(package_get_all_source_files SRCS PUBLIC_HEADERS PRIVATE_HEADERS)
  string(TOUPPER ${PROJECT_NAME} _project)

  unset(_tmp_srcs)
  unset(_tmp_public_headers)
  unset(_tmp_private_headers)

  package_get_all_activated_packages(_activated_list)
  foreach(_pkg_name ${_activated_list})
    _package_get_source_files(${_pkg_name}
      _pkg_srcs
      _pkg_public_headers
      _pkg_private_headers
      )
    list(APPEND _tmp_srcs ${_pkg_srcs})
    list(APPEND _tmp_public_headers ${_pkg_public_headers})
    list(APPEND _tmp_private_headers ${_pkg_private_headers})
  endforeach()

  set(${SRCS}            ${_tmp_srcs}            PARENT_SCOPE)
  set(${PUBLIC_HEADERS}  ${_tmp_public_headers}  PARENT_SCOPE)
  set(${PRIVATE_HEADERS} ${_tmp_private_headers} PARENT_SCOPE)
endfunction()


# ------------------------------------------------------------------------------
# Get include directories
# ------------------------------------------------------------------------------
function(package_get_all_include_directories inc_dirs)
  set(_tmp)
  package_get_all_activated_packages(_activated_list)

  foreach(_pkg_name ${_activated_list})
    foreach(_type SRCS PUBLIC_HEADERS PRIVATE_HEADERS)
      foreach(_file ${${_pkg_name}_${_type}})
        get_filename_component(_path "${_file}" PATH)
        list(APPEND _tmp "${_path}")
      endforeach()
    endforeach()
  endforeach()

  if(_tmp)
    list(REMOVE_DUPLICATES _tmp)
  endif()

  set(${inc_dirs} ${_tmp} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Get external libraries informations
# ------------------------------------------------------------------------------
function(package_get_all_external_informations)
  cmake_parse_arguments(_opt "" "PRIVATE_INCLUDE;INTERFACE_INCLUDE;LIBRARIES" "" ${ARGN})

  foreach(_type PRIVATE INTERFACE)
    if(_opt_${_type}_INCLUDE)
      _package_get_variable_for_external_dependencies(INCLUDE_DIR ${_type} tmp_INCLUDE_DIR)
      foreach(_dir ${tmp_INCLUDE_DIR})
        string(FIND "${_dir}" "${CMAKE_CURRENT_SOURCE_DIR}" _pos)
        if(NOT _pos EQUAL -1)
          list(REMOVE_ITEM tmp_INCLUDE_DIR ${_dir})
        endif()
      endforeach()

      set(${_opt_${_type}_INCLUDE} ${tmp_INCLUDE_DIR} PARENT_SCOPE)
    endif()
  endforeach()

  if(_opt_LIBRARIES)
    _package_get_variable_for_external_dependencies(LIBRARIES PRIVATE tmp_LIBRARIES)
    _package_get_variable_for_external_dependencies(LIBRARIES INTERFACE tmp_LIBRARIES_INTERFACE)
    set(${_opt_LIBRARIES} ${tmp_LIBRARIES} ${tmp_LIBRARIES_INTERFACE} PARENT_SCOPE)
  endif()
endfunction()

# ------------------------------------------------------------------------------
# Get export list for all activated packages
# ------------------------------------------------------------------------------
function(package_get_all_export_list export_list)
  _package_get_variable_for_activated(EXPORT_LIST _tmp)
  set(${export_list} ${_tmp} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Get definitions like external projects
# ------------------------------------------------------------------------------
function(package_get_all_definitions definitions)
  _package_get_variable_for_activated(OPTION_NAME _tmp)
  set(${definitions} ${_tmp} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Get extra dependencies like external projects
# ------------------------------------------------------------------------------
function(package_get_all_extra_dependencies deps)
  _package_get_variable_for_activated(EXTRA_DEPENDENCY _tmp)
  set(${deps} ${_tmp} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Get extra infos
# ------------------------------------------------------------------------------
function(package_get_all_test_folders TEST_DIRS)
  _package_get_variable_for_activated(TEST_FOLDER _tmp)
  set(${TEST_DIRS} ${_tmp} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Get compilation flags
# ------------------------------------------------------------------------------
function(package_get_all_compilation_flags LANG FLAGS)
  _package_get_variable_for_activated(COMPILE_${LANG}_FLAGS _tmp_flags)
  string(REPLACE ";" " " _flags "${_tmp_flags}")
  set(${FLAGS} ${_flags} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Documentation informations
# ------------------------------------------------------------------------------
function(package_get_all_documentation_files doc_files)
  set(_tmp_DOC_FILES)

  package_get_all_activated_packages(_activated_list)
  foreach(_pkg_name ${_activated_list})
    _package_get_manual_folder(${_pkg_name} _doc_dir)
    _package_get_documentation_files(${_pkg_name} _doc_files)

    foreach(_doc_file ${_doc_files})
      list(APPEND _tmp_DOC_FILES ${_doc_dir}/${_doc_file})
    endforeach()
  endforeach()

  if(_tmp_DOC_FILES)
    list(REMOVE_DUPLICATES _tmp_DOC_FILES)
  endif()

  set(${doc_files} ${_tmp_DOC_FILES} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Get package systems dependencies
# ------------------------------------------------------------------------------
function(package_get_all_package_system_dependency system deps)
  string(TOUPPER ${system} _u_system)
  _package_get_variable_for_activated(PACKAGE_SYSTEM_${_u_system} _tmp)
  set(${deps} ${_tmp} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# List packages
# ------------------------------------------------------------------------------
function(package_get_all_activated_packages activated_list)
  package_get_project_variable(ACTIVATED_PACKAGE_LIST _activated_list)
  set(${activated_list} ${_activated_list} PARENT_SCOPE)
endfunction()

function(package_get_all_deactivated_packages deactivated_list)
  package_get_project_variable(DEACTIVATED_PACKAGE_LIST _deactivated_list)
  set(${deactivated_list} ${_deactivated_list} PARENT_SCOPE)
endfunction()

function(package_get_all_packages packages_list)
  package_get_project_variable(ALL_PACKAGES_LIST _packages_list)
  set(${packages_list} ${_packages_list} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# List all the needed features
# ------------------------------------------------------------------------------
function(package_get_all_features_public features)
  _package_get_variable_for_activated(FEATURES_PUBLIC _tmp)
  set(${features} ${_tmp} PARENT_SCOPE)
endfunction()

function(package_get_all_features_private features)
  _package_get_variable_for_activated(FEATURES_PRIVATE _tmp)
  set(${features} ${_tmp} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Callbacks
# ------------------------------------------------------------------------------
function(package_on_enabled_script pkg script)
  package_get_name(${pkg} _pkg_name)
  _package_on_enable_script(${_pkg_name} "${script}")
endfunction()

# ------------------------------------------------------------------------------
# list all the packages in the PACKAGE_FOLDER
# extra packages can be given with an EXTRA_PACKAGE_FOLDER
# <package_folder>/<package>.cmake
#
# Extra packages folder structure
# <extra_package_folder>/<package>/package.cmake
#                                 /src
#                                 /test
#                                 /manual
#
# ------------------------------------------------------------------------------
function(package_list_packages PACKAGE_FOLDER)
  cmake_parse_arguments(_opt_pkg
    "NO_AUTO_COMPILE_FLAGS"
    "SOURCE_FOLDER;EXTRA_PACKAGES_FOLDER;TEST_FOLDER;MANUAL_FOLDER"
    ""
    ${ARGN})

  string(TOUPPER ${PROJECT_NAME} _project)

  # Cleaning some states to start correctly
  package_get_all_packages(_already_loaded_pkg)
  foreach(_pkg_name ${_already_loaded_pkg})
    _package_unset_extra_dependencies(${_pkg_name})
    _package_unset_dependencies(${_pkg_name} PRIVATE)
    _package_unset_dependencies(${_pkg_name} INTERFACE)
    _package_unset_activated(${_pkg_name})
  endforeach()


  if(_opt_pkg_SOURCE_FOLDER)
    set(_src_folder "${_opt_pkg_SOURCE_FOLDER}")
  else()
    set(_src_folder "src/")
  endif()

  get_filename_component(_abs_src_folder ${_src_folder} ABSOLUTE)

  if(_opt_pkg_TEST_FOLDER)
    set(_test_folder "${_opt_pkg_TEST_FOLDER}")
  else()
    set(_test_folder "test/")
  endif()

  if(_opt_pkg_MANUAL_FOLDER)
    set(_manual_folder "${_opt_pkg_MANUAL_FOLDER}")
  else()
    set(_manual_folder "doc/manual")
  endif()

  if(_opt_pkg_NO_AUTO_COMPILE_FLAGS)
    package_set_project_variable(NO_AUTO_COMPILE_FLAGS TRUE)
  else()
    package_set_project_variable(NO_AUTO_COMPILE_FLAGS FALSE)
  endif()


  get_filename_component(_abs_test_folder ${_test_folder} ABSOLUTE)
  get_filename_component(_abs_manual_folder ${_manual_folder} ABSOLUTE)

  # check all the packages in the <package_folder>
  file(GLOB _package_list "${PACKAGE_FOLDER}/*.cmake")

  set(_package_files)
  foreach(_pkg ${_package_list})
    get_filename_component(_basename ${_pkg} NAME)
    if(NOT _basename MATCHES "^\\.#.*")
      list(APPEND _package_files ${_basename})
    endif()
  endforeach()

  if(_package_files)
    list(SORT _package_files)
  endif()

  # check all packages
  set(_packages_list_all)
  foreach(_pkg_file ${_package_files})
    string(REGEX REPLACE "[0-9]+_" "" _pkg_file_stripped ${_pkg_file})
    string(REGEX REPLACE "\\.cmake" "" _pkg ${_pkg_file_stripped})

    set(_current_src_folder "${_abs_src_folder}" CACHE INTERNAL "" FORCE)
    set(_current_test_folder "${_abs_test_folder}" CACHE INTERNAL "" FORCE)
    set(_current_manual_folder "${_abs_manual_folder}" CACHE INTERNAL "" FORCE)

    include("${PACKAGE_FOLDER}/${_pkg_file}")

    unset(_current_src_folder CACHE)
    unset(_current_test_folder CACHE)
    unset(_current_manual_folder CACHE)
  endforeach()

  # check the extra_packages if they exists
  if(_opt_pkg_EXTRA_PACKAGES_FOLDER)
    file(GLOB _extra_package_list RELATIVE
      "${_opt_pkg_EXTRA_PACKAGES_FOLDER}" "${_opt_pkg_EXTRA_PACKAGES_FOLDER}/*")
    foreach(_pkg ${_extra_package_list})
      if(EXISTS "${_opt_pkg_EXTRA_PACKAGES_FOLDER}/${_pkg}/package.cmake")

        package_get_name(${_pkg} _pkg_name)

        _package_set_filename(${_pkg_name}
          "${_opt_pkg_EXTRA_PACKAGES_FOLDER}/${_pkg}/package.cmake")

        set(_current_src_folder "${_opt_pkg_EXTRA_PACKAGES_FOLDER}/${_pkg}/src" CACHE INTERNAL "" FORCE)

        if(EXISTS "${_opt_pkg_EXTRA_PACKAGES_FOLDER}/${_pkg}/test")
          set(_current_test_folder "${_opt_pkg_EXTRA_PACKAGES_FOLDER}/${_pkg}/test" CACHE INTERNAL "" FORCE)
        endif()

        if(EXISTS "${_opt_pkg_EXTRA_PACKAGES_FOLDER}/${_pkg}/manual")
          set(_current_manual_folder "${_opt_pkg_EXTRA_PACKAGES_FOLDER}/${_pkg}/manual" CACHE INTERNAL "" FORCE)
        endif()

        list(APPEND _extra_pkg_src_folders "${_opt_pkg_EXTRA_PACKAGES_FOLDER}/${_pkg}/src")

        include("${_opt_pkg_EXTRA_PACKAGES_FOLDER}/${_pkg}/package.cmake")

        unset(_current_src_folder CACHE)
        unset(_current_test_folder CACHE)
        unset(_current_manual_folder CACHE)
      endif()
    endforeach()
  endif()

  _package_build_rdependencies()
  _package_load_packages()
  _package_check_files_exists()
  _package_check_files_registered(${_abs_src_folder} ${_extra_pkg_src_folders})

  # Load boost components if boost was loaded
  package_is_activated(Boost _ret)
  if(_ret)
    _package_load_boost_components()
  endif()
endfunction()

# ------------------------------------------------------------------------------
# macro to include internal/external packages packages
# package_declare(<package real name>
#                 [EXTERNAL] [META] [ADVANCED] [NOT_OPTIONAL]
#                 [DESCRIPTION <description>] [DEFAULT <default_value>]
#                 [DEPENDS <pkg> ...]
#                 [BOOST_COMPONENTS <pkg> ...]
#                 [EXTRA_PACKAGE_OPTIONS <opt> ...]
#                 [COMPILE_FLAGS <lang> <flags>]
#                 [SYSTEM <bool> [ <script_to_compile> ]]
#                 [FEATURES_PUBLIC <feature> ...]
#                 [FEATURES_PRIVATE <feature> ...])
# ------------------------------------------------------------------------------
function(package_declare pkg)
  package_get_name(${pkg} _pkg_name)
  _package_set_real_name(${_pkg_name} ${pkg})
  _package_set_filename(${_pkg_name} "${CMAKE_CURRENT_LIST_FILE}")

  _package_set_sources_folder(${_pkg_name} "${_current_src_folder}")

  _package_variable_unset(SRCS ${_pkg_name})
  _package_variable_unset(PUBLIC_HEADERS ${_pkg_name})
  _package_variable_unset(PRIVATE_HEADERS ${_pkg_name})

  if(_current_test_folder)
    _package_set_tests_folder(${_pkg_name} "${_current_test_folder}")
  endif()

  if(_current_manual_folder)
    _package_set_manual_folder(${_pkg_name} "${_current_manual_folder}")
  endif()

  package_get_project_variable(ALL_PACKAGES_LIST _tmp_pkg_list)
  list(APPEND _tmp_pkg_list ${_pkg_name})
  list(REMOVE_DUPLICATES _tmp_pkg_list)
  package_set_project_variable(ALL_PACKAGES_LIST ${_tmp_pkg_list})

  set(_options
    EXTERNAL
    NOT_OPTIONAL
    META
    ADVANCED
    EXCLUDE_FROM_ALL)
  set(_one_valued_options
    DEFAULT
    DESCRIPTION)
  set(_multi_valued_options
    DEPENDS
    EXTRA_PACKAGE_OPTIONS
    COMPILE_FLAGS
    BOOST_COMPONENTS
    SYSTEM
    FEATURES_PUBLIC
    FEATURES_PRIVATE)

  cmake_parse_arguments(_opt_pkg
    "${_options}"
    "${_one_valued_options}"
    "${_multi_valued_options}"
    ${ARGN})

  if(_opt_pkg_UNPARSED_ARGUMENTS)
    message("You gave to many arguments while registering the package ${pkg} \"${_opt_pkg_UNPARSED_ARGUMENTS}\"")
  endif()

  # set the nature
  if(_opt_pkg_EXTERNAL)
    _package_set_nature(${_pkg_name} "external")
  elseif(_opt_pkg_META)
    _package_set_nature(${_pkg_name} "meta")
  else()
    _package_set_nature(${_pkg_name} "internal")
  endif()

  _package_declare_option(${_pkg_name})

  # set description
  if(_opt_pkg_DESCRIPTION)
    _package_set_description(${_pkg_name} ${_opt_pkg_DESCRIPTION})
  else()
    _package_set_description(${_pkg_name} "")
  endif()

  _package_get_option_name(${_pkg_name} _option_name)
  _package_get_description(${_pkg_name} _description)

  # get the default value
  if(DEFINED _opt_pkg_DEFAULT)
    set(_default ${_opt_pkg_DEFAULT})
  else()
    if(_opt_pkg_NOT_OPTIONAL)
      set(_default ON)
    else()
      set(_default OFF)
    endif()
  endif()

  # set the option if needed
  if(_opt_pkg_NOT_OPTIONAL)
    _package_get_nature(${_pkg_name} _nature)
    _package_set_nature(${_pkg_name} "${_nature}_not_optional")
    set(${_option_name} ${_default} CACHE INTERNAL "${_description}" FORCE)
    mark_as_advanced(${_option_name})
  else()
    option(${_option_name} "${_description}" ${_default})
    if(_opt_pkg_ADVANCED OR _opt_pkg_EXTERNAL)
      mark_as_advanced(${_option_name})
    endif()
  endif()

  # Set the option for third-partie that can be compiled as an ExternalProject
  if(DEFINED _opt_pkg_SYSTEM)
    list(LENGTH _opt_pkg_SYSTEM _length)
    list(GET _opt_pkg_SYSTEM 0 _bool)
    _package_set_system_option(${_pkg_name} ${_bool})
    if(_length GREATER 1)
      list(GET _opt_pkg_SYSTEM 1 _script)
      _package_set_system_script(${_pkg_name} ${_script})
    endif()
  endif()

  # set the dependecies
  if(_opt_pkg_DEPENDS)
    set(_deps_types PRIVATE PUBLIC INTERFACE)
    cmake_parse_arguments(_pkg_deps
      ""
      ""
      "${_deps_types}"
      ${_opt_pkg_DEPENDS})

    list(APPEND _pkg_deps_PRIVATE ${_pkg_deps_UNPARSED_ARGUMENTS})
    foreach(_type ${_deps_types})
      set(_depends)
      foreach(_dep ${_pkg_deps_${_type}})
	package_get_name(${_dep} _dep_pkg_name)
	list(APPEND _depends ${_dep_pkg_name})
      endforeach()
      _package_add_dependencies(${_pkg_name} ${_type} ${_depends})
    endforeach()
  endif()

  # keep the extra option for the future find package
  if(_opt_pkg_EXTRA_PACKAGE_OPTIONS)
    _package_set_find_package_extra_options(${_pkg_name} "${_opt_pkg_EXTRA_PACKAGE_OPTIONS}")
  endif()

  # register the compilation flags
  if(_opt_pkg_COMPILE_FLAGS)
    set(_languages C CXX Fortran)
    cmake_parse_arguments(_compile_flags
      "" "" "${_languages}"
      ${_opt_pkg_COMPILE_FLAGS}
      )


    # this is done to maintain backward compatibility
    if(_compile_flags_UNPARSED_ARGUMENTS)
      set(_compile_flags_CXX ${_compile_flags_UNPARSED_ARGUMENTS})
    endif()

    foreach(_lang ${_languages})
      if(_compile_flags_${_lang})
        _package_set_compile_flags(${_pkg_name} ${_lang} ${_compile_flags_${_lang}})
      else()
        _package_unset_compile_flags(${_pkg_name} ${_lang})
      endif()
    endforeach()
  endif()

  # set the boost dependencies
  if(_opt_pkg_BOOST_COMPONENTS)
    _package_set_boost_component_needed(${_pkg_name} "${_opt_pkg_BOOST_COMPONENTS}")
  endif()

  set(_variables FEATURES_PUBLIC FEATURES_PRIVATE EXCLUDE_FROM_ALL)
  foreach(_variable ${_variables})
    if(_opt_pkg_${_variable})
      _package_set_variable(${_variable} ${_pkg_name} "${_opt_pkg_${_variable}}")
    endif()
  endforeach()
endfunction()

# ------------------------------------------------------------------------------
# declare the source files of a given package
#
# package_declare_sources(<package> <list of sources>
#                         SOURCES <source file> ...
#                         PUBLIC_HEADER <header file> ...
#                         PRIVATE_HEADER <header file> ...)
# ------------------------------------------------------------------------------
function(package_declare_sources pkg)
  package_get_name(${pkg} _pkg_name)

  # get 3 lists, if none of the options given try to distinguish the different lists
  cmake_parse_arguments(_opt_pkg
    ""
    ""
    "SOURCES;PUBLIC_HEADERS;PRIVATE_HEADERS"
    ${ARGN})

  set(_tmp_srcs     ${_opt_pkg_SOURCES})
  set(_tmp_pub_hdrs ${_opt_pkg_PUBLIC_HEADER})
  set(_tmp_pri_hdrs ${_opt_pkg_PRIVATE_HEADERS})

  foreach(_file ${_opt_pkg_UNPARSED_ARGUMENTS})
    if(${_file} MATCHES ".*inline.*\\.cc")
      list(APPEND _tmp_pub_hdrs ${_file})
    elseif(${_file} MATCHES ".*\\.h+")
      list(APPEND _tmp_pub_hdrs ${_file})
    else()
      list(APPEND _tmp_srcs ${_file})
    endif()
  endforeach()

  _package_get_sources_folder(${_pkg_name} _src_folder)

  foreach(_type _srcs _pub_hdrs _pri_hdrs)
    set(${_type})
    foreach(_file ${_tmp${_type}})
      # get the full name
      set(_full_path "${_src_folder}/${_file}")
      list(APPEND ${_type} "${_full_path}")
    endforeach()
  endforeach()

  set(${_pkg_name}_SRCS "${_srcs}"
    CACHE INTERNAL "List of sources files" FORCE)
  set(${_pkg_name}_PUBLIC_HEADERS  "${_pub_hdrs}"
    CACHE INTERNAL "List of public header files" FORCE)
  set(${_pkg_name}_PRIVATE_HEADERS "${_pri_hdrs}"
    CACHE INTERNAL "List of private header files" FORCE)
endfunction()

# ------------------------------------------------------------------------------
