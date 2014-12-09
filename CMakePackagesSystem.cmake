#===============================================================================
# @file   CMakePackagesSystem.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
#
# @date creation: Thu Dec 20 2012
# @date last modification: Wed Sep 10 2014
#
# @brief  Set of macros used by akantu to handle the package system
#
# @section LICENSE
#
# Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#===============================================================================

include(CMakeParseArguments)

#===============================================================================
# Package Management
#===============================================================================
if(__CMAKE_PACKAGES_SYSTEM)
  return()
endif()
set(__CMAKE_PACKAGES_SYSTEM TRUE)


include(CMakeDebugMessages)
cmake_register_debug_message_module(PackagesSystem)

#===============================================================================
option(AUTO_MOVE_UNKNOWN_FILES
  "Give to cmake the permission to move the unregistered files to the ${PROJECT_SOURCE_DIR}/tmp directory" FALSE)

mark_as_advanced(AUTO_MOVE_UNKNOWN_FILES)

# ==============================================================================
# Accessors
# ==============================================================================
# Package name
function(package_get_name pkg pkg_name)
  string(TOUPPER ${PROJECT_NAME} _project)
  string(REPLACE "-" "_" _str_pkg "${pkg}")
  string(TOUPPER ${_str_pkg} _u_package)
  set(${pkg_name} ${_project}_PKG_${_u_package} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Real name
function(package_get_real_name pkg_name real_name)
  set(${real_name} ${${pkg_name}} PARENT_SCOPE)
endfunction()

function(_package_set_real_name pkg_name real_name)
  set(${pkg_name} ${real_name} CACHE INTERNAL "" FORCE)
endfunction()

# ------------------------------------------------------------------------------
# Option name
function(package_get_option_name pkg_name opt_name)
  string(TOUPPER "${PROJECT_NAME}" _project)
  package_get_real_name(${pkg_name} _real_name)
  string(TOUPPER "${_real_name}" _u_package)

  package_get_nature(${pkg_name} _nature)

  if(${_nature} MATCHES "internal" OR ${_nature} MATCHES "meta")
    set(${opt_name} ${_project}_${_u_package} PARENT_SCOPE)
  elseif(${_nature} MATCHES "external")
    set(${opt_name} ${_project}_USE_${_u_package} PARENT_SCOPE)
  else()
    set(${opt_name} UNKNOWN_NATURE_${_project}_${_u_package} PARENT_SCOPE)
  endif()
endfunction()

# ------------------------------------------------------------------------------
# Set if system package or compile external lib
function(_package_set_system_option pkg_name default)
  string(TOUPPER "${PROJECT_NAME}" _project)
  package_get_real_name(${pkg_name} _real_name)
  string(TOUPPER "${_real_name}" _u_package)

  option(${_project}_USE_SYSTEM_${_u_package}
    "Should akantu compile the third-party: \"${_real_name}\"" ${default})
  mark_as_advanced(${_project}_USE_SYSTEM_${_u_package})
endfunction()

function(package_use_system pkg_name use)
  string(TOUPPER "${PROJECT_NAME}" _project)
  package_get_real_name(${pkg_name} _real_name)
  string(TOUPPER "${_real_name}" _u_package)
  if(DEFINED ${_project}_USE_SYSTEM_${_u_package})
    set(${use} ${${_project}_USE_SYSTEM_${_u_package}} PARENT_SCOPE)
  else()
    set(${use} TRUE)
  endif()
endfunction()

# ------------------------------------------------------------------------------
# Nature
function(_package_set_nature pkg_name NATURE)
  set(${pkg_name}_NATURE ${NATURE} CACHE INTERNAL "" FORCE)
endfunction()

function(package_get_nature pkg_name NATURE)
  if(${pkg_name}_NATURE)
    set(${NATURE} ${${pkg_name}_NATURE} PARENT_SCOPE)
  else()
    set(${NATURE} "unknown" PARENT_SCOPE)
  endif()
endfunction()

# ------------------------------------------------------------------------------
# Description
function(_package_set_description pkg_name DESC)
  set(${pkg_name}_DESC ${DESC} CACHE INTERNAL "" FORCE)
endfunction()

function(package_get_description pkg_name DESC)
  if(${pkg_name}_DESC)
    set(${DESC} ${${pkg_name}_DESC} PARENT_SCOPE)
  else()
    message("No description set for the package ${${pkg_name}}")
  endif()
endfunction()

# ------------------------------------------------------------------------------
# Package file name
function(_package_set_filename pkg_name FILE)
  set(${pkg_name}_FILE ${FILE} CACHE INTERNAL "" FORCE)
endfunction()

function(package_get_filename pkg_name FILE)
  if(${pkg_name}_FILE)
    set(${FILE} ${${pkg_name}_FILE} PARENT_SCOPE)
  else()
    message(ERROR "No filename set for the package ${${pkg_name}}")
  endif()
endfunction()

# ------------------------------------------------------------------------------
# Source folder
function(_package_set_sources_folder pkg_name src_folder)
  set(${pkg_name}_SRCS_FOLDER ${src_folder} CACHE INTERNAL "" FORCE)
endfunction()

function(package_get_sources_folder pkg_name src_folder)
  set(${src_folder} ${${pkg_name}_SRCS_FOLDER} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Test folder
function(_package_set_tests_folder pkg_name test_folder)
  set(${pkg_name}_TESTS_FOLDER ${test_folder} CACHE INTERNAL "" FORCE)
endfunction()

function(package_get_tests_folder pkg_name test_folder)
  set(${test_folder} ${${pkg_name}_TESTS_FOLDER} PARENT_SCOPE)
endfunction()


# ------------------------------------------------------------------------------
# Extra option for the find_package
function(_package_set_extra_options pkg_name)
  set(${pkg_name}_OPTIONS ${ARGN}
    CACHE INTERNAL "Extra option for the fin_package function" FORCE)
endfunction()

function(package_get_extra_options pkg_name options)
  set(${options} "${${pkg_name}_OPTIONS}" PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Compilation flags
function(_package_set_compile_flags pkg_name)
  set(${pkg_name}_COMPILE_FLAGS ${ARGN}
    CACHE INTERNAL "Additional compile flags" FORCE)
endfunction()

function(package_get_compile_flags pkg_name flags)
  set(${flags} "${${pkg_name}_COMPILE_FLAGS}" PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Include dir
function(package_set_include_dir pkg_name)
  package_get_real_name(${pkg_name} _real_name)
  set(${pkg_name}_INCLUDE_DIR "${ARGN}"
    CACHE INTERNAL "Include folder for the package ${_real_name}" FORCE)
endfunction()

function(package_get_include_dir pkg_name include_dir)
  set(${include_dir} ${${pkg_name}_INCLUDE_DIR} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Libraries
function(package_set_libraries pkg_name)
  package_get_real_name(${pkg_name} _real_name)
  set(${pkg_name}_LIBRARIES "${ARGN}"
    CACHE INTERNAL "Libraries for the package ${_real_name}" FORCE)
endfunction()

function(package_get_libraries pkg_name libraries)
  set(${libraries} ${${pkg_name}_LIBRARIES} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Extra dependencies like custom commands of ExternalProject
function(package_add_extra_dependency PACKAGE)
  package_get_name(${PACKAGE} _pkg_name)
  set(_tmp_dep ${${pkg_name}_EXTRA_DEPENDS})
  list(APPEND _tmp_dep ${ARGN})
  list(REMOVE_DUPLICATES _tmp_dep)
  set(${pkg_name}_EXTRA_DEPENDS ${_tmp_dep} CACHE INTERNAL "External dependencies" FORCE)
endfunction()

function(package_get_extra_dependency PACKAGE deps)
  package_get_name(${PACKAGE} _pkg_name)
  set(${deps} ${${pkg_name}_EXTRA_DEPENDS} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Activate/deactivate
function(package_activate pkg_name)
  set(${pkg_name}_STATE ON CACHE INTERNAL "" FORCE)
endfunction()

function(package_deactivate pkg_name)
  set(${pkg_name}_STATE OFF CACHE INTERNAL "" FORCE)
endfunction()

function(package_is_activated pkg_name _act)
  if(${pkg_name}_STATE)
    set(${_act} TRUE PARENT_SCOPE)
  else()
    set(${_act} FALSE PARENT_SCOPE)
  endif()
endfunction()

# ------------------------------------------------------------------------------
# Direct dependencies
function(_package_set_dependencies pkg_name)
  set(${pkg_name}_DEPENDENCIES "${ARGN}"
    CACHE INTERNAL "List of dependencies for package ${_opt_name}" FORCE)
endfunction()

function(package_get_dependencies pkg_name dependencies)
  set(${dependencies} "${${pkg_name}_DEPENDENCIES}" PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Functions to handle reverse dependencies
# ------------------------------------------------------------------------------
function(_package_set_rdependencies pkg_name)
  set(${pkg_name}_RDEPENDENCIES "${ARGN}"
    CACHE INTERNAL "Dependencies ON with package ${pkg_name}" FORCE)
endfunction()

function(package_get_rdependencies pkg_name RDEPENDENCIES)
  set(${RDEPENDENCIES} "${${pkg_name}_RDEPENDENCIES}" PARENT_SCOPE)
endfunction()

function(_package_add_rdependency pkg_name rdep)
  # store the reverse dependency
  set(_rdeps ${${pkg_name}_RDEPENDENCIES})
  list(APPEND _rdeps ${rdep})
  list(REMOVE_DUPLICATES _rdeps)
  _package_set_rdependencies(${pkg_name} ${_rdeps})
endfunction()

function(_package_remove_rdependency pkg_name rdep)
  set(_rdeps ${${pkg_name}_RDEPENDENCIES})
  list(FIND _rdeps ${rdep} pos)
  if(NOT pos EQUAL -1)
    list(REMOVE_AT _rdeps ${pos})
    _package_set_rdependencies(${pkg_name} ${_rdeps})
  endif()
endfunction()

# ------------------------------------------------------------------------------
# Function to handle forcing dependencies (Package turn ON that enforce their
# dependencies ON)
# ------------------------------------------------------------------------------
function(_package_set_fdependencies pkg_name)
  set(${pkg_name}_FDEPENDENCIES "${ARGN}"
    CACHE INTERNAL "Dependencies ON with package ${pkg_name}" FORCE)
endfunction()

function(package_get_fdependencies pkg_name fdependencies)
  set(${fdependencies} "${${pkg_name}_FDEPENDENCIES}" PARENT_SCOPE)
endfunction()

function(_package_add_fdependency pkg_name fdep)
  # store the enforcing dependency
  set(_fdeps ${${pkg_name}_FDEPENDENCIES})
  list(APPEND _fdeps ${fdep})
  list(REMOVE_DUPLICATES _fdeps)
  _package_set_fdependencies(${pkg_name} ${_fdeps})
endfunction()

function(_package_remove_fdependency pkg_name fdep)
  set(_fdeps ${${pkg_name}_FDEPENDENCIES})
  list(FIND _fdeps ${fdep} pos)
  if(NOT pos EQUAL -1)
    list(REMOVE_AT _fdeps ${pos})
    _package_set_fdependencies(${pkg_name} ${_fdeps})
  endif()
endfunction()

# ==============================================================================
# Internal functions
# ==============================================================================

# ------------------------------------------------------------------------------
# Build the reverse dependencies from the dependencies
# ------------------------------------------------------------------------------
function(_package_build_rdependencies)
  string(TOUPPER ${PROJECT_NAME} _project)

  # set empty lists
  foreach(_pkg_name ${${_project}_ALL_PACKAGES_LIST})
    set(${_pkg_name}_rdeps)
  endforeach()

  # fill the dependencies list
  foreach(_pkg_name ${${_project}_ALL_PACKAGES_LIST})
    package_get_dependencies(${_pkg_name} _deps)
    foreach(_dep_name ${_deps})
      list(APPEND ${_dep_name}_rdeps ${_pkg_name})
    endforeach()
  endforeach()

  # clean and set the reverse dependencies
  foreach(_pkg_name ${${_project}_ALL_PACKAGES_LIST})
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
  string(TOUPPER ${PROJECT_NAME} _project)

  # Activate the dependencies of activated package and generate an ordered list
  # of dependencies
  set(ordered_loading_list)
  foreach(_pkg_name ${${_project}_ALL_PACKAGES_LIST})
    _package_load_dependencies_package(${_pkg_name} ordered_loading_list)
  endforeach()

  # Load the packages in the propoer order
  foreach(_pkg_name ${ordered_loading_list})
    package_use_system(${_pkg_name} _use_system)
    if(_use_system)
      _package_load_package(${_pkg_name})
    endif()
  endforeach()

  # generates the activated and unactivated lists of packages
  set(_packages_activated)
  set(_packages_deactivated)
  foreach(_pkg_name ${${_project}_ALL_PACKAGES_LIST})
    package_is_activated(${_pkg_name} _act)
    if(_act)
      list(APPEND _packages_activated ${_pkg_name})
    else()
      list(APPEND _packages_deactivated ${_pkg_name})
    endif()
  endforeach()

  # generate the list usable by the calling code
  set(${_project}_ACTIVATED_PACKAGE_LIST "${_packages_activated}"
    CACHE INTERNAL "List of activated packages" FORCE)
  set(${_project}_DEACTIVATED_PACKAGE_LIST "${_packages_deactivated}"
    CACHE INTERNAL "List of deactivated packages" FORCE)
endfunction()



# ------------------------------------------------------------------------------
# This load an external package and recursively all its dependencies
# ------------------------------------------------------------------------------
function(_package_load_dependencies_package pkg_name loading_list)
  # Get packages informations
  package_get_option_name(${pkg_name} _pkg_option_name)
  package_get_dependencies(${pkg_name} _dependencies)

  # handle the dependencies
  foreach(_dep_name ${_dependencies})
    package_get_description(${_dep_name} _dep_desc)
    package_get_option_name(${_dep_name} _dep_option_name)

    package_get_fdependencies(${_dep_name} _fdeps)
    if(${_pkg_option_name})
      if("${_fdeps}" STREQUAL "")
        set(${_dep_name}_OLD ${${_dep_option_name}} CACHE INTERNAL "" FORCE)
      endif()

      # set the option to on
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
  package_get_compile_flags(${pkg_name} _pkg_comile_flags)

  # if package option is on add it in the list
  if(${_pkg_option_name})
    list(FIND ${loading_list} ${pkg_name} _pos)
    if(_pos EQUAL -1)
      set(_tmp_loading_list ${${loading_list}})
      list(APPEND _tmp_loading_list ${pkg_name})
      set(${loading_list} "${_tmp_loading_list}" PARENT_SCOPE)
    endif()

    #add the comilation flags if needed
    if(_pkg_comile_flags)
      add_flags(cxx ${_pkg_comile_flags})
    endif()
  else()
    package_deactivate(${pkg_name})

    #remove the comilation flags if needed
    if(_pkg_comile_flags)
      remove_flags(cxx ${_pkg_comile_flags})
    endif()
  endif()
endfunction()

# ------------------------------------------------------------------------------
# Load the package if it is an external one
# ------------------------------------------------------------------------------
function(_package_load_package pkg_name)
  package_get_option_name(${pkg_name} _pkg_option_name)

  # load the package if it is an external
  package_get_nature(${pkg_name} _nature)
  if(${_nature} MATCHES "external")
    package_get_extra_options(${pkg_name} _options)
    if(_options)
      cmake_parse_arguments(_opt_pkg "" "LANGUAGE" "PREFIX;FOUND;ARGS" ${_options})
      if(_opt_pkg_UNPARSED_ARGUMENTS)
	message("You passed too many options for the find_package related to ${${pkg_name}} \"${_opt_pkg_UNPARSED_ARGUMENTS}\"")
      endif()
    endif()

    if(_opt_pkg_LANGUAGE)
      foreach(_language ${_opt_pkg_LANGUAGE})
        enable_language(${_language})
      endforeach()
    endif()

    package_get_real_name(${pkg_name} _real_name)

    # find the package
    find_package(${_real_name} REQUIRED ${_opt_pkg_ARGS})

    # check if the package is found
    if(_opt_pkg_PREFIX)
      set(_package_prefix ${_opt_pkg_PREFIX})
    else()
      string(TOUPPER ${${pkg_name}} _u_package)
      set(_package_prefix ${_u_package})
    endif()

    foreach(_prefix ${_package_prefix})
      if(${_prefix}_FOUND OR _opt_pkg_FOUND)
	# Add the in the definition list
        list(APPEND ${_project}_DEFINITIONS ${_option_name})

	# Generate the include dir for the package
        if(DEFINED ${_prefix}_INCLUDE_DIRS)
          package_set_include_dir(${_pkg_name} ${${_prefix}_INCLUDE_DIRS})
        elseif(DEFINED ${_prefix}_INCLUDE_DIR)
	  package_set_include_dir(${_pkg_name} ${${_prefix}_INCLUDE_DIR})
        elseif(DEFINED ${_prefix}_INCLUDE_PATH)
	  package_set_include_dir(${_pkg_name} ${${_prefix}_INCLUDE_PATH})
        endif()

	# Generate the libraries for the package
	package_set_libraries(${_pkg_name} ${${_prefix}_LIBRARIES})

        package_activate(${_pkg_name})
      else()
	package_deactivate(${_pkg_name})
      endif()
    endforeach()
  elseif(${_nature} MATCHES "internal")
    package_activate(${_pkg_name})
  endif()
endfunction()

# ------------------------------------------------------------------------------
# Sanity check functions
# ------------------------------------------------------------------------------
function(_package_check_files_exists)
  string(TOUPPER ${PROJECT_NAME} _project)

  set(_message FALSE)

  foreach(_pkg_name ${${_project}_ALL_PACKAGES_LIST})
    set(_pkg_files
      ${${_pkg_name}_SRCS}
      ${${_pkg_name}_PUBLIC_HEADERS}
      ${${_pkg_name}_PRIVATE_HEADERS}
      )

    package_get_real_name(${_pkg_name} _real_name)

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
  # generates a file list of registered files
  foreach(_pkg_name ${${_project}_ALL_PACKAGES_LIST})
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
  # check only sources files ine the source folders
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
    file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/tmp/)
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
        file(RENAME ${_file} ${PROJECT_SOURCE_DIR}/tmp/${_file_name})
      endif()

      file(APPEND ${PROJECT_BINARY_DIR}/missing_files_in_packages "${_file}
")
    endforeach()

    if(AUTO_MOVE_UNKNOWN_FILES)
      message(SEND_ERROR "The files where moved in the followinf folder ${PROJECT_SOURCE_DIR}/tmp/")
    endif()
    message(SEND_ERROR "Please register them in the good package or clean the sources")
  endif()
endfunction()

# ==============================================================================
# User Functions
# ==============================================================================

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
    ""
    "SOURCE_FOLDER;EXTRA_PACKAGES_FOLDER;TEST_FOLDER;MANUAL_FOLDER"
    ""
    ${ARGN})

  string(TOUPPER ${PROJECT_NAME} _project)

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

  get_filename_component(_abs_test_folder ${_test_folder} ABSOLUTE)

  # check all the packages in the <package_folder>
  file(GLOB _package_list "${PACKAGE_FOLDER}/*.cmake")

  set(_package_files)
  foreach(_pkg ${_package_list})
    get_filename_component(_basename ${_pkg} NAME)
    list(APPEND _package_files ${_basename})
  endforeach()

  if(_package_files)
    list(SORT _package_files)
  endif()

  # check all packages
  set(_packages_list_all)
  foreach(_pkg_file ${_package_files})
    string(REGEX REPLACE "[0-9]+_" "" _pkg_file_stripped ${_pkg_file})
    string(REGEX REPLACE "\\.cmake" "" _pkg ${_pkg_file_stripped})

    package_get_name(${_pkg} _pkg_name)
    _package_set_filename(${_pkg_name} "${PACKAGE_FOLDER}/${_pkg_file}")
    _package_set_sources_folder(${_pkg_name} "${_abs_src_folder}")

    _package_set_tests_folder(${_pkg_name}
      "${_abs_test_folder}")

    list(APPEND _packages_list_all ${_pkg_name})
    include("${PACKAGE_FOLDER}/${_pkg_file}")
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

	_package_set_sources_folder(${_pkg_name}
	  "${_opt_pkg_EXTRA_PACKAGES_FOLDER}/${_pkg}/src")

	if(EXISTS "${_opt_pkg_EXTRA_PACKAGES_FOLDER}/${_pkg}/test")
	  _package_set_tests_folder(${_pkg_name}
	    "${_opt_pkg_EXTRA_PACKAGES_FOLDER}/${_pkg}/test")
	endif()

	list(APPEND _extra_pkg_src_folders "${_opt_pkg_EXTRA_PACKAGES_FOLDER}/${_pkg}/src")

	list(APPEND _packages_list_all ${_pkg_name})
	include("${_opt_pkg_EXTRA_PACKAGES_FOLDER}/${_pkg}/package.cmake")
      endif()
    endforeach()
  endif()

  # Store the list of packages
  string(TOUPPER ${PROJECT_NAME} _project)
  set(${_project}_ALL_PACKAGES_LIST ${_packages_list_all}
    CACHE INTERNAL "List of available packages" FORCE)

  _package_build_rdependencies()
  _package_load_packages()
  _package_check_files_exists()
  _package_check_files_registered(${_abs_src_folder} ${_extra_pkg_src_folders})
endfunction()

# ------------------------------------------------------------------------------
# macro to include internal/external packages packages
# package_declare(<package real name>
#                 [EXTERNAL] [META] [ADVANCED] [NOT_OPTIONAL]
#                 [DESCRIPTION <description>] [DEFAULT <default_value>]
#                 [DEPENDS <pkg> ...]
#                 [EXTRA_PACKAGE_OPTIONS <opt> ...]
#                 [COMPILE_FLAGS <flags>])
# ------------------------------------------------------------------------------
function(package_declare PACKAGE)
  package_get_name(${PACKAGE} _pkg_name)
  _package_set_real_name(${_pkg_name} ${PACKAGE})

  cmake_parse_arguments(_opt_pkg
    "EXTERNAL;NOT_OPTIONAL;META;ADVANCED"
    "DEFAULT;DESCRIPTION;SYSTEM"
    "DEPENDS;EXTRA_PACKAGE_OPTIONS;COMPILE_FLAGS"
    ${ARGN})

  if(_opt_pkg_UNPARSED_ARGUMENTS)
    message("You gave to many arguments while registering the package ${PACKAGE} \"${_opt_pkg_UNPARSED_ARGUMENTS}\"")
  endif()

  # set description
  if(_opt_pkg_DESCRIPTION)
    _package_set_description(${_pkg_name} ${_opt_pkg_DESCRIPTION})
  else()
    _package_set_description(${_pkg_name} "")
  endif()

  # set the nature
  if(_opt_pkg_EXTERNAL)
    _package_set_nature(${_pkg_name} "external")
  elseif(_opt_pkg_META)
    _package_set_nature(${_pkg_name} "meta")
  else()
    _package_set_nature(${_pkg_name} "internal")
  endif()

  package_get_option_name(${_pkg_name} _option_name)
  package_get_description(${_pkg_name} _description)

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
    package_get_nature(${_pkg_name} _nature)
    _package_set_nature(${_pkg_name} "${_nature}_not_optional")
    set(${_option_name} ${_default} CACHE INTERNAL "${_description}" FORCE)
  else()
    option(${_option_name} "${_description}" ${_default})
    if(_opt_pkg_ADVANCED OR _opt_pkg_EXTERNAL)
      mark_as_advanced(${_option_name})
    endif()
  endif()

  # Set the option for third-partie that can be compiled as an ExternalProject
  if(DEFINED _opt_pkg_SYSTEM)
    _package_set_system_option(${_pkg_name} ${_opt_pkg_SYSTEM})
  endif()

  # set the dependecies
  if(_opt_pkg_DEPENDS)
    set(_depends)
    foreach(_dep ${_opt_pkg_DEPENDS})
      package_get_name(${_dep} _dep_pkg_name)
      list(APPEND _depends ${_dep_pkg_name})
    endforeach()
    _package_set_dependencies(${_pkg_name} ${_depends})
    package_get_dependencies(${_pkg_name} _deps)
  endif()

  # keep the extra option for the future find package
  if(_opt_pkg_EXTRA_PACKAGE_OPTIONS)
    _package_set_extra_options(${_pkg_name} "${_opt_pkg_EXTRA_PACKAGE_OPTIONS}")
  endif()

  # register the compilation flags
  if(_opt_pkg_COMPILE_FLAGS)
    _package_set_compile_flags(${_pkg_name} "${_opt_pkg_COMPILE_FLAGS}")
  endif()

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

  # get 3 lists, if non of the options given try to distinguish the different lists
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

  package_get_sources_folder(${_pkg_name} _src_folder)

  foreach(_type _srcs _pub_hdrs _pri_hdrs)
    set(${_type})
    foreach(_file ${_tmp${_type}})
      # get the full name
      list(APPEND ${_type} "${_src_folder}/${_file}")
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
# get the list of source files
# ------------------------------------------------------------------------------
function(package_get_source_files SRCS PUBLIC_HEADERS PRIVATE_HEADERS)
  string(TOUPPER ${PROJECT_NAME} _project)

  set(tmp_SRCS)
  set(tmp_PUBLIC_HEADERS)
  set(tmp_PRIVATE_HEADERS)

  foreach(_pkg_name ${${_project}_ACTIVATED_PACKAGE_LIST})
    foreach(_type SRCS PUBLIC_HEADERS PRIVATE_HEADERS)
      foreach(_file ${${_pkg_name}_${_type}})
	string(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/" "" _rel_file "${_file}")
	list(APPEND tmp_${_type} "${_rel_file}")
      endforeach()
    endforeach()
  endforeach()

  set(${SRCS}            ${tmp_SRCS}            PARENT_SCOPE)
  set(${PUBLIC_HEADERS}  ${tmp_PUBLIC_HEADERS}  PARENT_SCOPE)
  set(${PRIVATE_HEADERS} ${tmp_PRIVATE_HEADERS} PARENT_SCOPE)
endfunction()


# ------------------------------------------------------------------------------
# Get include directories
# ------------------------------------------------------------------------------
function(package_get_include_directories inc_dirs)
  string(TOUPPER ${PROJECT_NAME} _project)

  set(_tmp)
  foreach(_pkg_name ${${_project}_ACTIVATED_PACKAGE_LIST})
    foreach(_type SRCS PUBLIC_HEADERS PRIVATE_HEADERS)
      foreach(_file ${${_pkg_name}_${_type}})
	get_filename_component(_path "${_file}" PATH)
	list(APPEND _tmp "${_path}")
      endforeach()
    endforeach()
  endforeach()

  list(REMOVE_DUPLICATES _tmp)

  set(${inc_dirs} ${_tmp} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Get external libraries informations
# ------------------------------------------------------------------------------
function(package_get_external_informations INCLUDE_DIR LIBRARIES)
  set(tmp_INCLUDE_DIR)
  set(tmp_LIBRARIES)
  foreach(_pkg_name ${${_project}_ACTIVATED_PACKAGE_LIST})
    package_get_nature(${_pkg_name} _nature)
    if(${_nature} MATCHES "external")
      package_get_include_dir(${_pkg_name} _inc)
      package_get_libraries  (${_pkg_name} _lib)

      list(APPEND tmp_INCLUDE_DIR ${_inc})
      list(APPEND tmp_LIBRARIES   ${_lib})
    endif()
  endforeach()

  set(${INCLUDE_DIR} ${tmp_INCLUDE_DIR} PARENT_SCOPE)
  set(${LIBRARIES}   ${tmp_LIBRARIES}   PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Get extra dependencies like external projects
# ------------------------------------------------------------------------------
function(package_get_all_extra_dependency DEPS)
  set(_tmp_DEPS)
  foreach(_pkg_name ${${_project}_ACTIVATED_PACKAGE_LIST})
    package_get_extra_dependency(${_pkg_name} _dep)
    list(APPEND _tmp_DEPS ${_dep})
  endforeach()
  if(_tmp_DEPS)
    list(REMOVE_DUPLICATES _tmp_DEPS)
  endif()
  set(${DEPS} ${_tmp_DEPS} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Get extra dependencies like external projects
# ------------------------------------------------------------------------------
function(package_get_all_test_folders TEST_DIRS)
  set(_tmp_TEST_DIRS)
  foreach(_pkg_name ${${_project}_ACTIVATED_PACKAGE_LIST})
    package_get_tests_folder(${_pkg_name} _test_dir)
    list(APPEND _tmp_TEST_DIRS ${_test_dir})
  endforeach()
  if(_tmp_TEST_DIRS)
    list(REMOVE_DUPLICATES _tmp_TEST_DIRS)
  endif()
  set(${TEST_DIRS} ${_tmp_TEST_DIRS} PARENT_SCOPE)
endfunction()
