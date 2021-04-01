#===============================================================================
# @file   CMakePackagesSystemPrivateFunctions.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Sat Jul 18 2015
# @date last modification: Wed Jan 20 2016
#
# @brief  Set of macros used by the package system, internal functions
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

if (DEFINED CMAKE_PACKAGES_SYSTEM_PRIVATE_FUNCTIONS_LOADED)
  return()
endif()
set(CMAKE_PACKAGES_SYSTEM_PRIVATE_FUNCTIONS_LOADED TRUE)

# ==============================================================================
# "Private" Accessors
# ==============================================================================
# ------------------------------------------------------------------------------
# Real name
# ------------------------------------------------------------------------------
function(_package_get_real_name pkg_name real_name)
  set(${real_name} ${${pkg_name}} PARENT_SCOPE)
endfunction()

function(_package_set_real_name pkg_name real_name)
  set(${pkg_name} ${real_name} CACHE INTERNAL "" FORCE)
endfunction()

# ------------------------------------------------------------------------------
# Option name
# ------------------------------------------------------------------------------
function(_package_declare_option pkg_name)
  string(TOUPPER "${PROJECT_NAME}" _project)
  _package_get_real_name(${pkg_name} _real_name)
  string(TOUPPER "${_real_name}" _u_package)

  _package_get_nature(${pkg_name} _nature)

  if(${_nature} MATCHES "internal" OR ${_nature} MATCHES "meta")
    set(_opt_name ${_project}_${_u_package})
  elseif(${_nature} MATCHES "external")
    set(_opt_name ${_project}_USE_${_u_package})
  else()
    set(_opt_name UNKNOWN_NATURE_${_project}_${_u_package})
  endif()

  _package_set_variable(OPTION_NAME ${pkg_name} ${_opt_name})
endfunction()

function(_package_get_option_name pkg_name opt_name)
  _package_get_variable(OPTION_NAME ${pkg_name} _opt_name)
  set(${opt_name} ${_opt_name} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Set if system package or compile external lib
# ------------------------------------------------------------------------------
function(_package_set_system_option pkg_name default)
  string(TOUPPER "${PROJECT_NAME}" _project)
  _package_get_real_name(${pkg_name} _real_name)
  string(TOUPPER "${_real_name}" _u_package)

  set(${_project}_USE_SYSTEM_${_u_package} ${default} CACHE STRING
    "Should akantu compile the third-party: \"${_real_name}\"")
  mark_as_advanced(${_project}_USE_SYSTEM_${_u_package})
  set_property(CACHE ${_project}_USE_SYSTEM_${_u_package} PROPERTY STRINGS ON OFF AUTO)
endfunction()

function(_package_use_system pkg_name use)
  string(TOUPPER "${PROJECT_NAME}" _project)
  _package_get_real_name(${pkg_name} _real_name)
  string(TOUPPER "${_real_name}" _u_package)
  if(DEFINED ${_project}_USE_SYSTEM_${_u_package})
    if(${${_project}_USE_SYSTEM_${_u_package}} MATCHES "(ON|AUTO)")
      set(${use} TRUE PARENT_SCOPE)
    else()
      set(${use} FALSE PARENT_SCOPE)
    endif()
  else()
    set(${use} TRUE PARENT_SCOPE)
  endif()
endfunction()

function(_package_has_system_fallback pkg_name fallback)
  string(TOUPPER "${PROJECT_NAME}" _project)
  _package_get_real_name(${pkg_name} _real_name)
  string(TOUPPER "${_real_name}" _u_package)

  if(DEFINED ${_project}_USE_SYSTEM_${_u_package})
    if(${${_project}_USE_SYSTEM_${_u_package}} MATCHES "AUTO")
      set(${fallback} TRUE PARENT_SCOPE)
    else()
      set(${fallback} FALSE PARENT_SCOPE)
    endif()
  else()
    set(${fallback} FALSE PARENT_SCOPE)
  endif()
endfunction()

function(_package_set_system_script pkg_name script)
  _package_set_variable(COMPILE_SCRIPT ${pkg_name} "${script}")
endfunction()

function(_package_add_third_party_script_variable pkg_name var)
  _package_set_variable(VARIABLE_${var} ${pkg_name} "${ARGN}")
  set(${var} ${ARGN} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
function(_package_load_third_party_script pkg_name)
  if(${pkg_name}_COMPILE_SCRIPT)
    # set the stored variable
    get_cmake_property(_all_vars VARIABLES)
    foreach(_var ${_all_vars})
      if(_var MATCHES "^${pkg_name}_VARIABLE_.*")
        string(REPLACE "${pkg_name}_VARIABLE_" "" _orig_var "${_var}")
        set(${_orig_var} ${${_var}})
      endif()
    endforeach()

    _package_get_real_name(${pkg_name} _name)
    string(TOUPPER "${_name}" _u_name)

    _package_get_option_name(${pkg_name} _opt_name)
    if(${_opt_name}_VERSION)
      set(_version "${${_opt_name}_VERSION}")
      set(${_u_name}_VERSION "${_version}" CACHE INTERNAL "" FORCE)
    elseif(${_u_name}_VERSION)
      set(_version "${${_u_name}_VERSION}")
    endif()

    # load the script
    include(ExternalProject)
    include(${${pkg_name}_COMPILE_SCRIPT})

    if(${_u_name}_LIBRARIES)
      _package_set_libraries(${pkg_name} ${${_u_name}_LIBRARIES})
      list(APPEND _required_vars ${_u_name}_LIBRARIES)
    endif()
    
    if(${_u_name}_INCLUDE_DIR)
      _package_set_include_dir(${pkg_name} ${${_u_name}_INCLUDE_DIR})
      list(APPEND _required_vars ${_u_name}_INCLUDE_DIR)
    endif()

    include(FindPackageHandleStandardArgs)
    if (NOT _required_vars)
      message(FATAL_ERROR "The package ${_name} does not define any of the variables ${_u_name}_INCLUDE_DIR nor ${_u_name}_LIBRARIES")
    endif()

    if(CMAKE_VERSION VERSION_GREATER 2.8.12)
      find_package_handle_standard_args(${_name}
        REQUIRED_VARS ${_required_vars}
        VERSION_VAR _version
        FAIL_MESSAGE "Something was not configured by a the third-party script for ${_name}"
        )
    else()
      find_package_handle_standard_args(${_name}
        "Something was not configured by a the third-party script for ${_name}"
        ${_required_vars}
        )
    endif()
  endif()
  set(${pkg_name}_USE_SYSTEM_PREVIOUS FALSE CACHE INTERNAL "" FORCE)
endfunction()

# ------------------------------------------------------------------------------
# Nature
# ------------------------------------------------------------------------------
function(_package_set_nature pkg_name nature)
  _package_set_variable(NATURE ${pkg_name} ${nature})
endfunction()

function(_package_get_nature pkg_name nature)
  _package_get_variable(NATURE ${pkg_name} _nature "unknown")
  set(${nature} ${_nature} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Description
# ------------------------------------------------------------------------------
function(_package_set_description pkg_name desc)
  _package_set_variable(DESC ${pkg_name} ${desc})
endfunction()

function(_package_get_description pkg_name desc)
  _package_get_variable(DESC ${pkg_name} _desc "No description set for the package ${${pkg_name}} (${pkg_name})")
  set(${desc} ${_desc} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Package file name
# ------------------------------------------------------------------------------
function(_package_set_filename pkg_name file)
  _package_set_variable(FILE ${pkg_name} ${file})
endfunction()

function(_package_get_filename pkg_name file)
  _package_get_variable(FILE ${pkg_name} _file "No filename set for the package ${${pkg_name}}")
  set(${file} ${_file} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Source folder
# ------------------------------------------------------------------------------
function(_package_set_sources_folder pkg_name src_folder)
  _package_set_variable(SRC_FOLDER ${pkg_name} ${src_folder})
endfunction()

function(_package_get_sources_folder pkg_name src_folder)
  _package_get_variable(SRC_FOLDER ${pkg_name} _src_folder)
  set(${src_folder} ${_src_folder} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Test folder
# ------------------------------------------------------------------------------
function(_package_set_tests_folder pkg_name test_folder)
  _package_set_variable(TEST_FOLDER ${pkg_name} ${test_folder})
endfunction()

function(_package_get_tests_folder pkg_name test_folder)
  _package_get_variable(TEST_FOLDER ${pkg_name} _test_folder)
  set(${test_folder} ${_test_folder} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Manual folder
# ------------------------------------------------------------------------------
function(_package_set_manual_folder pkg_name manual_folder)
  _package_set_variable(MANUAL_FOLDER ${pkg_name} ${manual_folder})
endfunction()

function(_package_get_manual_folder pkg_name manual_folder)
  _package_get_variable(MANUAL_FOLDER ${pkg_name} _manual_folder)
  set(${manual_folder} ${_manual_folder} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Extra option for the find_package
# ------------------------------------------------------------------------------
function(_package_set_find_package_extra_options pkg_name)
  _package_set_variable(FIND_PKG_OPTIONS ${pkg_name} ${ARGN})
endfunction()

function(_package_get_find_package_extra_options pkg_name options)
  _package_get_variable(FIND_PKG_OPTIONS ${pkg_name} _options)
  set(${options} ${_options} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Compilation flags
# ------------------------------------------------------------------------------
function(_package_set_compile_flags pkg_name lang)
  _package_set_variable(COMPILE_${lang}_FLAGS ${pkg_name} ${ARGN})
endfunction()

function(_package_unset_compile_flags pkg_name lang)
  _package_variable_unset(COMPILE_${lang}_FLAGS ${pkg_name})
endfunction()

function(_package_get_compile_flags pkg_name lang flags)
  _package_get_variable(COMPILE_${lang}_FLAGS ${pkg_name} _tmp_flags)
  string(REPLACE ";" " " _flags "${_tmp_flags}")
  set(${flags} "${_flags}" PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Include dir
# ------------------------------------------------------------------------------
function(_package_set_include_dir pkg_name)
  _package_set_variable(INCLUDE_DIR ${pkg_name} ${ARGN})
endfunction()

function(_package_get_include_dir pkg_name include_dir)
  _package_get_variable(INCLUDE_DIR ${pkg_name} _include_dir "")
  set(${include_dir} ${_include_dir} PARENT_SCOPE)
endfunction()

function(_package_add_include_dir pkg_name)
  _package_add_to_variable(INCLUDE_DIR ${pkg_name} ${ARGN})
endfunction()

# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------
function(_package_set_libraries pkg_name)
  _package_set_variable(LIBRARIES ${pkg_name} ${ARGN})
endfunction()

function(_package_get_libraries pkg_name libraries)
  _package_get_variable(LIBRARIES ${pkg_name} _libraries "")
  set(${libraries} ${_libraries} PARENT_SCOPE)
endfunction()

function(_package_add_libraries pkg_name)
  _package_add_to_variable(LIBRARIES ${pkg_name} ${ARGN})
endfunction()

# ------------------------------------------------------------------------------
# Extra dependencies like custom commands of ExternalProject
# ------------------------------------------------------------------------------
function(_package_add_extra_dependency pkg_name)
  _package_add_to_variable(EXTRA_DEPENDENCY ${pkg_name} ${ARGN})
endfunction()

function(_package_rm_extra_dependency pkg_name dep)
  _package_remove_from_variable(EXTRA_DEPENDENCY ${pkg_name} ${dep})
endfunction()

function(_package_set_extra_dependencies pkg)
  _package_set_variable(EXTRA_DEPENDENCY ${pkg_name} ${ARGN})
endfunction()

function(_package_get_extra_dependencies pkg deps)
  _package_get_variable(EXTRA_DEPENDENCY ${pkg_name} _deps "")
  set(${deps} ${_deps} PARENT_SCOPE)
endfunction()

function(_package_unset_extra_dependencies pkg_name)
  _package_variable_unset(EXTRA_DEPENDENCY ${pkg_name})
endfunction()

# ------------------------------------------------------------------------------
# Activate/deactivate
# ------------------------------------------------------------------------------
function(_package_activate pkg_name)
  _package_set_variable(STATE ${pkg_name} ON)
endfunction()

function(_package_deactivate pkg_name)
  _package_set_variable(STATE ${pkg_name} OFF)
endfunction()

function(_package_is_activated pkg_name act)
  _package_get_variable(STATE ${pkg_name} _state OFF)
  if(_state)
    set(${act} TRUE PARENT_SCOPE)
  else()
    set(${act} FALSE PARENT_SCOPE)
  endif()
endfunction()

function(_package_is_deactivated pkg_name act)
  _package_get_variable(STATE ${pkg_name} _state OFF)
  if(NOT _state)
    set(${act} TRUE PARENT_SCOPE)
  else()
    set(${act} FALSE PARENT_SCOPE)
  endif()
endfunction()

function(_package_unset_activated pkg_name)
  _package_variable_unset(STATE ${pkg_name})
endfunction()

# ------------------------------------------------------------------------------
# Callbacks
# ------------------------------------------------------------------------------
function(_package_on_enable_script pkg_name script)
  string(TOLOWER "${pkg_name}" _l_pkg_name)
  set(_output_file "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${_l_pkg_name}.cmake")
  file(WRITE "${_output_file}"
    "${script}")
  _package_set_variable(CALLBACK_SCRIPT ${pkg_name} "${_output_file}")
endfunction()

function(_package_get_callback_script pkg_name filename)
  _package_get_variable(CALLBACK_SCRIPT ${pkg_name} _filename NOTFOUND)
  set(${filename} ${_filename} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Export list
# ------------------------------------------------------------------------------
function(_package_add_to_export_list pkg_name)
  _package_add_to_variable(EXPORT_LIST ${pkg_name} ${ARGN})
endfunction()

function(_package_remove_from_export_list pkg_name)
  _package_remove_from_variable(EXPORT_LIST ${pkg_name} ${ARGN})
endfunction()

function(_package_get_export_list pkg_name export_list)
  _package_get_variable(EXPORT_LIST ${pkg_name} _export_list)
  set(${export_list} ${_export_list} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Direct dependencies
# ------------------------------------------------------------------------------
function(_package_add_dependencies pkg_name type)
  _package_add_to_variable(DEPENDENCIES_${type} ${pkg_name} ${ARGN})
endfunction()

function(_package_get_dependencies pkg_name type dependencies)
  _package_get_variable(DEPENDENCIES_${type} ${pkg_name} _dependencies)
  set(${dependencies} ${_dependencies} PARENT_SCOPE)
endfunction()

function(_package_unset_dependencies pkg_name type)
  _package_variable_unset(DEPENDENCIES_${type} ${pkg_name})
endfunction()

function(_package_remove_dependency pkg_name type dep)
  set(_deps ${${pkg_name}_DEPENDENCIES_${type}})

  _package_get_fdependencies(${dep} _fdeps)

  # check if this is the last reverse dependency
  list(LENGTH _fdeps len)
  list(FIND _fdeps ${pkg_name} pos)
  if((len EQUAL 1) AND (NOT pos EQUAL -1))
    _package_get_description(${dep} _dep_desc)
    _package_get_option_name(${dep} _dep_option_name)
    set(${_dep_option_name} ${${dep}_OLD} CACHE BOOL "${_dep_desc}" FORCE)
    unset(${dep}_OLD CACHE)
  endif()

  # remove the pkg_name form the reverse dependency
  _package_remove_fdependency(${dep} ${pkg_name})

  list(FIND _deps ${dep} pos)
  if(NOT pos EQUAL -1)
    list(REMOVE_AT _deps ${pos})
    _package_set_variable(DEPENDENCIES_${type} ${pkg_name} ${_deps})
  endif()
endfunction()

# ------------------------------------------------------------------------------
# Functions to handle reverse dependencies
# ------------------------------------------------------------------------------
function(_package_set_rdependencies pkg_name)
  _package_set_variable(RDEPENDENCIES ${pkg_name} ${ARGN})
endfunction()

function(_package_get_rdependencies pkg_name rdependencies)
  _package_get_variable(RDEPENDENCIES ${pkg_name} _rdependencies)
  set(${rdependencies} ${_rdependencies} PARENT_SCOPE)
endfunction()

function(_package_add_rdependency pkg_name rdep)
  # store the reverse dependency
  _package_add_to_variable(RDEPENDENCIES ${pkg_name} ${rdep})
endfunction()

function(_package_remove_rdependency pkg_name rdep)
  _package_remove_from_variable(RDEPENDENCIES ${pkg_name} ${rdep})
endfunction()

# ------------------------------------------------------------------------------
# Function to handle forcing dependencies (Package turn ON that enforce their
# dependencies ON)
# ------------------------------------------------------------------------------
function(_package_set_fdependencies pkg_name)
  _package_set_variable(FDEPENDENCIES ${pkg_name} ${ARGN})
endfunction()

function(_package_get_fdependencies pkg_name fdependencies)
  _package_get_variable(FDEPENDENCIES ${pkg_name} _fdependencies)
  set(${fdependencies} ${_fdependencies} PARENT_SCOPE)
endfunction()

function(_package_add_fdependency pkg_name fdep)
  # store the reverse dependency
  _package_add_to_variable(FDEPENDENCIES ${pkg_name} ${fdep})
endfunction()

function(_package_remove_fdependency pkg_name fdep)
  _package_remove_from_variable(FDEPENDENCIES ${pkg_name} ${fdep})
endfunction()

# ------------------------------------------------------------------------------
# Exteral package system as apt rpm dependencies
# ------------------------------------------------------------------------------
function(_package_set_package_system_dependency pkg system)
  string(TOUPPER "${_system}" _u_system)
  _package_set_variable(PACKAGE_SYSTEM_${_u_system} ${_pkg_name} ${ARGN})
endfunction()


function(_package_get_package_system_dependency pkg system var)
  string(TOUPPER "${_system}" _u_system)
  _package_get_variable(PACKAGE_SYSTEM_${_u_system} ${_pkg_name} ${_deps})
  set(${var} ${_deps} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Documentation related functions
# ------------------------------------------------------------------------------
function(_package_set_documentation_files pkg_name)
  _package_set_variable(DOCUMENTATION_FILES ${pkg_name} ${ARGN})
endfunction()

function(_package_get_documentation_files pkg_name doc_files)
  _package_get_variable(DOCUMENTATION_FILES ${pkg_name} _doc_files "")
  set(${doc_files} ${_doc_files} PARENT_SCOPE)
endfunction()

function(_package_set_documentation pkg_name)
  # \n replaced by && and \\ by ££ to avoid cache problems
  set(_doc_str "")
  foreach(_str ${ARGN})
    set(_doc_str "${_doc_str}&&${_str}")
  endforeach()

  string(REPLACE "\\" "££" _doc_escaped "${_doc_str}")

  _package_set_variable(DOCUMENTATION ${pkg_name} "${_doc_str}")
endfunction()

function(_package_get_documentation pkg_name _doc)
  # \n replaced by && and \\ by ££ to avoid cache problems
  _package_get_variable(DOCUMENTATION ${pkg_name} _doc_tmp "")

  string(REPLACE "££" "\\" _doc_escaped "${_doc_tmp}")
  string(REPLACE "&&" "\n" _doc_newlines "${_doc_escaped}")
  set(${_doc} "${_doc_newlines}" PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Special boost thingy
# ------------------------------------------------------------------------------
function(_package_set_boost_component_needed pkg_name)
  _package_add_to_variable(BOOST_COMPONENTS_NEEDED ${pkg_name} ${ARGN})
  package_get_name(Boost _boost_pkg_name)
  _package_add_dependencies(${pkg_name} PRIVATE ${_boost_pkg_name})
endfunction()

function(_package_get_boost_component_needed pkg_name needed)
  _package_get_variable(BOOST_COMPONENTS_NEEDED ${pkg_name} _needed)
  set(${needed} ${_needed} PARENT_SCOPE)
endfunction()

function(_package_load_boost_components)
  string(TOUPPER ${PROJECT_NAME} _project)

  _package_get_variable_for_activated(BOOST_COMPONENTS_NEEDED _boost_needed_components)
  package_get_name(Boost _pkg_name)

  if(_boost_needed_components)
    message(STATUS "Looking for Boost liraries: ${_boost_needed_components}")
    foreach(_comp ${_boost_needed_components})
      find_package(Boost COMPONENTS ${_comp} QUIET)
      string(TOUPPER ${_comp} _u_comp)
      if(Boost_${_u_comp}_FOUND)
        message(STATUS "   ${_comp}: FOUND")
        package_set_project_variable(BOOST_${_u_comp} TRUE)

        # Generate the libraries for the package
        _package_add_libraries(${_pkg_name} ${Boost_${_u_comp}_LIBRARY})
      else()
        message(STATUS "   ${_comp}: NOT FOUND")
      endif()
    endforeach()
  endif()
endfunction()

# ------------------------------------------------------------------------------
# get the list of source files for a given package
# ------------------------------------------------------------------------------
function(_package_get_source_files pkg_name SRCS PUBLIC_HEADERS PRIVATE_HEADERS)
  string(TOUPPER ${PROJECT_NAME} _project)

  set(tmp_SRCS)
  set(tmp_PUBLIC_HEADERS)
  set(tmp_PRIVATE_HEADERS)

  foreach(_type SRCS PUBLIC_HEADERS PRIVATE_HEADERS)
    foreach(_file ${${pkg_name}_${_type}})
      string(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/" "" _rel_file "${_file}")
      list(APPEND tmp_${_type} "${_rel_file}")
    endforeach()
  endforeach()

  set(${SRCS}            ${tmp_SRCS}            PARENT_SCOPE)
  set(${PUBLIC_HEADERS}  ${tmp_PUBLIC_HEADERS}  PARENT_SCOPE)
  set(${PRIVATE_HEADERS} ${tmp_PRIVATE_HEADERS} PARENT_SCOPE)
endfunction()

# ==============================================================================
# Internal functions
# ==============================================================================
# ------------------------------------------------------------------------------
# Build the reverse dependencies from the dependencies
# ------------------------------------------------------------------------------
function(_package_build_rdependencies)
  package_get_all_packages(_pkg_list)

  # set empty lists
  foreach(_pkg_name ${_pkg_list})
    set(${_pkg_name}_rdeps)
  endforeach()

  # fill the dependencies list
  foreach(_pkg_name ${_pkg_list})
    _package_get_dependencies(${_pkg_name} PRIVATE _deps_priv)
    _package_get_dependencies(${_pkg_name} INTERFACE _deps_interface)
    foreach(_dep_name ${_deps_priv} ${_deps_interface})
      list(APPEND ${_dep_name}_rdeps ${_pkg_name})
    endforeach()
  endforeach()

  # clean and set the reverse dependencies
  foreach(_pkg_name ${_pkg_list})
    if(${_pkg_name}_rdeps)
      list(REMOVE_DUPLICATES ${_pkg_name}_rdeps)
      _package_set_rdependencies(${_pkg_name} ${${_pkg_name}_rdeps})
    endif()
  endforeach()
endfunction()

# ------------------------------------------------------------------------------
# This function resolve the dependance order run the needed find_packages
# ------------------------------------------------------------------------------
function(_package_load_packages)
  package_get_all_packages(_pkg_list)

  # Activate the dependencies of activated package and generate an ordered list
  # of dependencies
  set(ordered_loading_list)
  foreach(_pkg_name ${_pkg_list})
    _package_load_dependencies_package(${_pkg_name} ordered_loading_list)
  endforeach()

  # Load the packages in the propoer order
  foreach(_pkg_name ${ordered_loading_list})
    _package_get_option_name(${_pkg_name} _option_name)

    if(${_option_name})
      _package_load_package(${_pkg_name})
    else()
      # deactivate the packages than can already be deactivated
      _package_deactivate(${_pkg_name})
    endif()
  endforeach()

  # generates the activated and unactivated lists of packages
  set(_packages_activated)
  set(_packages_deactivated)
  foreach(_pkg_name ${_pkg_list})
    _package_is_activated(${_pkg_name} _active)
    set(_exclude FALSE)
    _package_get_variable(EXCLUDE_FROM_ALL ${_pkg_name} _exclude)
    if(_active AND NOT _exclude)
      list(APPEND _packages_activated ${_pkg_name})
    else()
      list(APPEND _packages_deactivated ${_pkg_name})
    endif()
  endforeach()

  # generate the list usable by the calling code
  package_set_project_variable(ACTIVATED_PACKAGE_LIST "${_packages_activated}")
  package_set_project_variable(DEACTIVATED_PACKAGE_LIST "${_packages_deactivated}")
endfunction()



# ------------------------------------------------------------------------------
# This load an external package and recursively all its dependencies
# ------------------------------------------------------------------------------
function(_package_load_dependencies_package pkg_name loading_list)
  # Get packages informations
  _package_get_option_name(${pkg_name} _pkg_option_name)
  _package_get_dependencies(${pkg_name} PRIVATE _deps_private)
  _package_get_dependencies(${pkg_name} INTERFACE _deps_interface)
  set(_dependencies ${_deps_private} ${_deps_interface})

  # handle the dependencies
  foreach(_dep_name ${_dependencies})
    _package_get_description(${_dep_name} _dep_desc)
    _package_get_option_name(${_dep_name} _dep_option_name)

    _package_get_fdependencies(${_dep_name} _fdeps)
    if(${_pkg_option_name})
      if("${_fdeps}" STREQUAL "")
        set(${_dep_name}_OLD ${${_dep_option_name}} CACHE INTERNAL "" FORCE)
      endif()

      # set the option to on
      set(_type BOOL)
      _package_get_nature(${_dep_name} _nature)
      if(_nature MATCHES ".*not_optional$")
	set(_type INTERNAL)
      endif()
      set(${_dep_option_name} ON CACHE BOOL "${_dep_desc}" FORCE)

      # store the reverse dependency
      _package_add_fdependency(${_dep_name} ${pkg_name})
    else()
      # check if this is the last reverse dependency
      list(LENGTH _fdeps len)
      list(FIND _fdeps ${pkg_name} pos)
      if((len EQUAL 1) AND (NOT pos EQUAL -1))
        set(${_dep_option_name} ${${_dep_name}_OLD} CACHE BOOL "${_dep_desc}" FORCE)
        unset(${_dep_name}_OLD CACHE)
      endif()

      # remove the pkg_name form the reverse dependency
      _package_remove_fdependency(${_dep_name} ${pkg_name})
    endif()

    # recusively load the dependencies
    _package_load_dependencies_package(${_dep_name} ${loading_list})
  endforeach()

  # get the compile flags
  _package_get_compile_flags(${pkg_name} CXX _pkg_compile_flags)
  package_get_project_variable(NO_AUTO_COMPILE_FLAGS _no_auto_compile_flags)

  # if package option is on add it in the list
  if(${_pkg_option_name})
    list(FIND ${loading_list} ${pkg_name} _pos)
    if(_pos EQUAL -1)
      set(_tmp_loading_list ${${loading_list}})
      list(APPEND _tmp_loading_list ${pkg_name})
      set(${loading_list} "${_tmp_loading_list}" PARENT_SCOPE)
    endif()

    #add the comilation flags if needed
    if(_pkg_compile_flags AND NOT _no_auto_compile_flags)
      add_flags(cxx ${_pkg_compile_flags})
    endif()
  else()
    #remove the comilation flags if needed
    if(_pkg_comile_flags AND NOT _no_auto_compile_flags)
      remove_flags(cxx ${_pkg_compile_flags})
    endif()
  endif()
endfunction()

# ------------------------------------------------------------------------------
# Load the package if it is an external one
# ------------------------------------------------------------------------------
function(_package_load_package pkg_name)
  # load the package if it is an external
  _package_get_nature(${pkg_name} _nature)
  if(${_nature} MATCHES "external")
    _package_use_system(${pkg_name} _use_system)

    set(_activated TRUE)
    if(_use_system)
      _package_load_external_package(${pkg_name} _activated)
    endif()

    _package_has_system_fallback(${pkg_name} _fallback)
    if((NOT _use_system) OR (_fallback AND (NOT _activated)))
      _package_load_third_party_script(${pkg_name})
      set(_activated TRUE)
    endif()

    if(_activated)
      _package_activate(${pkg_name})
    elseif()
      _package_deactivate(${pkg_name})
    endif()
  else()
    _package_activate(${pkg_name})
  endif()
endfunction()


# ------------------------------------------------------------------------------
# Load external packages
# ------------------------------------------------------------------------------
function(_package_load_external_package pkg_name activate)
  _package_get_real_name(${pkg_name} _real_name)
  string(TOUPPER ${_real_name} _u_package)

  if(NOT ${pkg_name}_USE_SYSTEM_PREVIOUS)
    #if system was off before clear the cache of preset variables
    get_cmake_property(_all_vars VARIABLES)
    foreach(_var ${_all_vars})
      if(_var MATCHES "^${_u_package}_.*")
        unset(${_var} CACHE)
      endif()
    endforeach()
    set(${pkg_name}_USE_SYSTEM_PREVIOUS TRUE CACHE INTERNAL "" FORCE)
  endif()

  _package_get_find_package_extra_options(${pkg_name} _options)
  if(_options)
    cmake_parse_arguments(_opt_pkg "" "LANGUAGE" "LINK_LIBRARIES;PREFIX;FOUND;ARGS" ${_options})
    if(_opt_pkg_UNPARSED_ARGUMENTS)
      message("You passed too many options for the find_package related to ${${pkg_name}} \"${_opt_pkg_UNPARSED_ARGUMENTS}\"")
    endif()
  endif()

  if(_opt_pkg_LANGUAGE)
    foreach(_language ${_opt_pkg_LANGUAGE})
      enable_language(${_language})
    endforeach()
  endif()

  # find the package
  set(_required REQUIRED)
  _package_has_system_fallback(${pkg_name} _fallback)
  if(_fallback)
    set(_required QUIET)
  endif()

  find_package(${_real_name} ${_required} ${_opt_pkg_ARGS})

  # check if the package is found
  if(_opt_pkg_PREFIX)
    set(_package_prefix ${_opt_pkg_PREFIX})
  else()
    set(_package_prefix ${_u_package})
  endif()

  set(_act FALSE)
  set(_prefix_to_consider)
  if(DEFINED _opt_pkg_FOUND)
    if(${_opt_pkg_FOUND})
      set(_act TRUE)
      set(_prefix_to_consider ${_package_prefix})
    endif()
  else()
    foreach(_prefix ${_package_prefix})
      if(${_prefix}_FOUND)
        set(_act TRUE)
        list(APPEND _prefix_to_consider ${_prefix})
      endif()
    endforeach()
  endif()

  if(_act)
    foreach(_prefix ${_prefix_to_consider})
      # Generate the include dir for the package
      if(DEFINED ${_prefix}_INCLUDE_DIRS)
        _package_set_include_dir(${pkg_name} ${${_prefix}_INCLUDE_DIRS})
      elseif(DEFINED ${_prefix}_INCLUDE_DIR)
        _package_set_include_dir(${pkg_name} ${${_prefix}_INCLUDE_DIR})
      elseif(DEFINED ${_prefix}_INCLUDE_PATH)
        _package_set_include_dir(${pkg_name} ${${_prefix}_INCLUDE_PATH})
      endif()

      # Generate the libraries for the package
      if(_opt_pkg_LINK_LIBRARIES)
        _package_set_libraries(${pkg_name} ${_opt_pkg_LINK_LIBRARIES})
      elseif(DEFINED ${_prefix}_LIBRARIES)
        _package_set_libraries(${pkg_name} ${${_prefix}_LIBRARIES})
      elseif(DEFINED ${_prefix}_LIBRARY)
        _package_set_libraries(${pkg_name} ${${_prefix}_LIBRARY})
      endif()
    endforeach()

    _package_get_callback_script(${pkg_name} _script_file)
    if(_script_file)
      include("${_script_file}")
    endif()
  endif()
  set(${activate} ${_act} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Sanity check functions
# ------------------------------------------------------------------------------
function(_package_check_files_exists)
  set(_message FALSE)

  package_get_all_packages(_pkg_list)
  foreach(_pkg_name ${_pkg_list})
    set(_pkg_files
      ${${_pkg_name}_SRCS}
      ${${_pkg_name}_PUBLIC_HEADERS}
      ${${_pkg_name}_PRIVATE_HEADERS}
      )

    _package_get_real_name(${_pkg_name} _real_name)

    foreach(_file ${_pkg_files})
      if(NOT EXISTS "${_file}")
        if(NOT _message)
          set(_message TRUE)
          message("This file(s) is(are) present in a package but are not present on disk.")
        endif()

        message(" PACKAGE ${_real_name} FILE ${_file}")
      endif()
    endforeach()
  endforeach()

  if(_message)
    message(SEND_ERROR "Please check the content of your packages to correct this warnings")
  endif()
endfunction()

# ------------------------------------------------------------------------------
function(_package_check_files_registered)
  set(_pkg_files)
  package_get_all_packages(_pkg_list)
  # generates a file list of registered files
  foreach(_pkg_name ${_pkg_list})
    list(APPEND _pkg_files
      ${${_pkg_name}_SRCS}
      ${${_pkg_name}_PUBLIC_HEADERS}
      ${${_pkg_name}_PRIVATE_HEADERS}
      )
  endforeach()

  # generates the list of files in the source folders
  set(_all_src_files)
  foreach(_src_folder ${ARGN})
    foreach(_ext "cc" "hh" "c" "h" "hpp")
      file(GLOB_RECURSE _src_files "${_src_folder}/*.${_ext}")
      list(APPEND _all_src_files ${_src_files})
    endforeach()
  endforeach()

  if(_all_src_files)
    list(REMOVE_DUPLICATES _all_src_files)
  endif()

  set(_not_registerd_files)
  # check only sources files in the source folders
  foreach(_src_folder ${ARGN})
    foreach(_file ${_all_src_files})
      if("${_file}" MATCHES  "${_src_folder}")
        list(FIND _pkg_files "${_file}" _index)
        if (_index EQUAL -1)
          list(APPEND _not_registerd_files ${_file})
        endif()
      endif()
    endforeach()
  endforeach()

  if(AUTO_MOVE_UNKNOWN_FILES)
    file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/unknown_files)
  endif()

  # warn the user and move the files if needed
  if(_not_registerd_files)
    if(EXISTS ${PROJECT_BINARY_DIR}/missing_files_in_packages)
      file(REMOVE ${PROJECT_BINARY_DIR}/missing_files_in_packages)
    endif()

    message("This files are present in the source folders but are not registered in any package")
    foreach(_file ${_not_registerd_files})
      message(" ${_file}")
      if(AUTO_MOVE_UNKNOWN_FILES)
        get_filename_component(_file_name ${_file} NAME)
        file(RENAME ${_file} ${PROJECT_BINARY_DIR}/unknown_files/${_file_name})
      endif()

      file(APPEND ${PROJECT_BINARY_DIR}/missing_files_in_packages "${_file}
")
    endforeach()

    if(AUTO_MOVE_UNKNOWN_FILES)
      message(SEND_ERROR "The files where moved in the followinf folder ${PROJECT_BINARY_DIR}/unknown_files\n
Please register them in the good package or clean the sources")
    else()
      message(SEND_ERROR "Please register them in the good package or clean the sources")
    endif()

  endif()
endfunction()

# ------------------------------------------------------------------------------
