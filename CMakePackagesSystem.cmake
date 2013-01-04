#===============================================================================
# @file   CMakePackagesSystem.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date   Tue Oct 16 14:05:02 2012
#
# @brief  Set of macros used by akantu to handle the package system
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
#===============================================================================

#===============================================================================
# Package Management
#===============================================================================
if(__CMAKE_PACKAGES_SYSTEM)
  return()
endif()
set(__CMAKE_PACKAGES_SYSTEM TRUE)

include(CMakeDebugMessages)
cmake_register_debug_message_module(PackagesSystem)

macro(package_pkg_name PKG PKG_NAME)
  string(TOUPPER ${PROJECT_NAME} _project)
  string(TOUPPER ${PKG} _u_package)
  set(${PKG_NAME} ${_project}_${_u_package})
endmacro()

macro(package_opt_name PKG OPT_NAME)
  string(TOUPPER ${PROJECT_NAME} _project)
  string(TOUPPER ${PKG} _u_package)
  set(${OPT_NAME} ${_project}_USE_${_u_package})
endmacro()

macro(package_desc_name PKG DESC_NAME)
  string(TOUPPER ${PROJECT_NAME} _project)
  string(TOUPPER ${PKG} _u_package)
  set(${DESC_NAME} ${_project}_DESC_${_u_package})
endmacro()

#===============================================================================
macro(add_all_packages package_dir src_dir)
  string(TOUPPER ${PROJECT_NAME} _project)
  cmake_debug_message(PackagesSystem "add_all_packages: PKG DIR : ${package_dir}")
  file(GLOB ${_project}_package_list "${package_dir}/*.cmake")

  set(_${_project}_src_dir ${src_dir})

  foreach(_pkg ${${_project}_package_list})
    get_filename_component(_basename ${_pkg} NAME)
    string(REGEX REPLACE "\\.cmake" "" _option_name ${_basename})
    string(TOUPPER "${_option_name}" _option_name)
    list(APPEND ${_project}_PACKAGE_SYSTEM_PACKAGES_NAMES_LIST_ALL ${_option_name})
  endforeach()

  cmake_debug_message(PackagesSystem "add_all_packages: PKG LIST : ${${_project}_PACKAGE_SYSTEM_PACKAGES_NAMES_LIST_ALL}")

  foreach(_pkg ${${_project}_PACKAGE_SYSTEM_PACKAGES_NAMES_LIST_ALL})
    cmake_debug_message(PackagesSystem "add_all_packages: including ${_pkg}")
    string(TOLOWER "${_pkg}" _pkg)
    include(${package_dir}/${_pkg}.cmake)
    package_pkg_name(${_pkg} _package_name)
    if (${_package_name})
      list(APPEND ${_project}_PACKAGE_SYSTEM_PACKAGES_ON ${_pkg})
    else (${_package_name})
      list(APPEND ${_project}_PACKAGE_SYSTEM_PACKAGES_OFF ${_pkg})
    endif()

    foreach(_file ${${_package_name}_FILES})
      list(APPEND ${_project}_release_all_files ${_file})
    endforeach()
  endforeach()

  cmake_debug_message(PackagesSystem "add_all_packages: ON  PKG : ${${_project}_PACKAGE_SYSTEM_PACKAGES_ON}")
  cmake_debug_message(PackagesSystem "add_all_packages: ALL RELEASE FILES LIST : ${_project}_release_all_files ${${_project}_release_all_files}")

  #check if there are some file in the release that are not registered in a package
  file(GLOB_RECURSE ${_project}_all_files "*.cc" "*.hh" "*.c" "*.h")
  cmake_debug_message(PackagesSystem "add_all_packages: ALL FILES LIST : ${_project}_all_files ${${_project}_all_files}")

  cmake_debug_message(PackagesSystem "add_all_packages: SOURCE DIR : ${_${_project}_src_dir}")
  foreach(_file ${${_project}_all_files})
    if("${_file}" MATCHES  "${_${_project}_src_dir}")
      file(RELATIVE_PATH __file "${_${_project}_src_dir}" ${_file})
      list(APPEND ${_project}_all_files_relatives ${__file})
    endif()
  endforeach()

  foreach(_file ${${_project}_all_files_relatives})
    if(NOT ${_file} MATCHES "test.*" AND NOT ${_file} MATCHES "third-party" AND NOT ${_file} MATCHES "doc")
      list(FIND ${_project}_release_all_files ${_file} _index)
      if (_index EQUAL -1)
	list(APPEND ${_project}_missing_files_in_packages ${_file})
        message("The file ${_file} is not registered in any package.")
        message("Please append the file in one of the files within directory ${PROJECT_SOURCE_DIR}/packages")
      endif()
    endif()
  endforeach()

  if(${_project}_missing_files_in_packages)
    message("A complete list of files missing in the packeges description can be found here: ${PROJECT_BINARY_DIR}/missing_files_in_packages")
    if(EXISTS ${PROJECT_BINARY_DIR}/missing_files_in_packages)
      file(REMOVE ${PROJECT_BINARY_DIR}/missing_files_in_packages)
    endif()
    foreach(_file ${_missing_files_in_packages})
      file(APPEND ${PROJECT_BINARY_DIR}/missing_files_in_packages "${_file}
")
    endforeach()
  endif()

  #check if there are some file in the package list that are not on the current directory
  foreach(_file ${${_project}_release_all_files})
    list(FIND ${_project}_all_files_relatives ${_file} _index)
    if (_index EQUAL -1)
      message("The file ${_file} is registered in packages but is not present in the source directory.")
    endif()
  endforeach()

  #construct list of files for unactivated packages
  foreach(_pkg ${${_project}_PACKAGE_SYSTEM_PACKAGES_OFF})
    package_pkg_name(${_pkg} _pkg_name)
    foreach(_file ${${_pkg_name}_FILES})
      cmake_debug_message(PackagesSystem "add_all_packages: ${_file}")
      list(APPEND ${_project}_EXCLUDE_SOURCE_FILE ${_file})
    endforeach()
  endforeach()

  #check dependencies
  foreach(_pkg ${${_project}_PACKAGE_SYSTEM_PACKAGES_OFF})
    # differentiate the file types
    cmake_debug_message(PackagesSystem "add_all_packages: DEPENDS PKG : ${_pkg}")
    cmake_debug_message(PackagesSystem "add_all_packages: DEPENDS LST : ${${_pkg}_DEPENDS}")
    package_pkg_name(${_pkg} _pkg_name)
    if (NOT "${${_pkg_name}_DEB_DEPEND}" STREQUAL "")
      set(deb_depend "${deb_depend}, ${${_pkg}_DEB_DEPEND}")
    endif()
  endforeach()
  set(${_project}_PACKAGE_SYSTEM_DEBIAN_PACKAGE_DEPENDS "${deb_depend}")
endmacro()

#===============================================================================
macro(generate_source_list_from_packages source_dir source_files inline_files headers_files include_dirs)
  string(TOUPPER ${PROJECT_NAME} _project)

  cmake_debug_message(PackagesSystem "generate_source_list_from_packages: SRC DIR : ${source_dir}")
  foreach(_pkg ${${_project}_PACKAGE_SYSTEM_PACKAGES_ON})
    # differentiate the file types
    package_pkg_name(${_pkg} _package_name)
    cmake_debug_message(PackagesSystem "generate_source_list_from_packages: PKG ${_package_name} FILES : ${${_package_name}_FILES}")
    foreach(_file ${${_package_name}_FILES})
      if(${_file} MATCHES ".*inline.*\\.cc")
	list(APPEND ${_package_name}_inlines ${_file})
      elseif(${_file} MATCHES ".*\\.h+")
        list(APPEND ${_package_name}_headers ${_file})
      else()
        list(APPEND ${_package_name}_srcs ${_file})
      endif()
    endforeach()

    # generates the include directory variable
    foreach(_file ${${_package_name}_headers})
      get_filename_component(_absolute_name ${_file} ABSOLUTE)
      get_filename_component(_include_dir ${_absolute_name} PATH)
      list(APPEND ${_package_name}_include_dirs ${_include_dir})
      list(REMOVE_DUPLICATES ${_package_name}_include_dirs)
    endforeach()

    # generate global lists for akantu to know what to build
    list(APPEND ${source_files} ${${_package_name}_srcs})
    list(APPEND ${inline_files} ${${_package_name}_inlines})
    list(APPEND ${headers_files} ${${_package_name}_headers})
    list(APPEND ${include_dirs}  ${${_package_name}_include_dirs})

    cmake_debug_message(PackagesSystem "generate_source_list_from_packages: PKG ${_package_name} SRCS : ${${_package_name}_srcs}")
    cmake_debug_message(PackagesSystem "generate_source_list_from_packages: PKG ${_package_name} INLINES : ${${_package_name}_inlines}")
    cmake_debug_message(PackagesSystem "generate_source_list_from_packages: PKG ${_package_name} HRDS : ${${_package_name}_headers}")
    cmake_debug_message(PackagesSystem "generate_source_list_from_packages: PKG ${_package_name} INCS : ${${_package_name}_include_dirs}")
  endforeach()

  cmake_debug_message(PackagesSystem "generate_source_list_from_packages: SRCS : ${${source_files}}")
  cmake_debug_message(PackagesSystem "generate_source_list_from_packages: HRDS : ${${headers_files}}")
  cmake_debug_message(PackagesSystem "generate_source_list_from_packages: INCS : ${${include_dirs}}")
endmacro()

#===============================================================================
# macro to include optional packages
macro(add_optional_external_package PACKAGE DESC DEFAULT)
  package_opt_name (${PACKAGE} _option_name)
  package_desc_name(${PACKAGE} _desc_name)
  set(${_desc_name} ${DESC})
  option(${_option_name} ${DESC} ${DEFAULT})
  _add_external_package(${PACKAGE} ${ARGN})
endmacro()

macro(add_external_package PACKAGE)
  package_opt_name (${PACKAGE} _option_name)
  set(${_option_name} ON)
  _add_external_package(${PACKAGE} ${ARGN})
endmacro()


macro(_add_external_package PACKAGE)
  string(TOUPPER ${PROJECT_NAME} _project)
  cmake_parse_arguments(_opt_pkg "" "LANGUAGE" "DEPENDS;PREFIX;FOUND;ARGS" ${ARGN})

  package_pkg_name (${PACKAGE} _pkg_name)
  package_opt_name (${PACKAGE} _option_name)

  cmake_debug_message(PackagesSystem "add_optional_package: Registering ${PACKAGE} ${DESC} -> ${_option_name}")

  if(_opt_pkg_PREFIX)
    set(_package_prefix ${_opt_pkg_PREFIX})
  else()
    string(TOUPPER ${PACKAGE} _u_package)
    set(_package_prefix ${_u_package})
  endif()

  if(${_option_name})
    if(_opt_pkg_LANGUAGE)
      foreach(_language ${_opt_pkg_LANGUAGE})
	cmake_debug_message(PackagesSystem "add_optional_package: Package ${PACKAGE} asked for language ${_language}")
	enable_language(${_language})
      endforeach()
    endif()

    foreach(_dep ${_opt_pkg_DEPENDS})
      add_external_package_dependencies(${PACKAGE} ${_dep})
    endforeach()

    find_package(${PACKAGE} REQUIRED ${_opt_pkg_ARGS})

    foreach(_prefix ${_package_prefix})
      if(${_prefix}_FOUND OR _opt_pkg_FOUND)
	list(APPEND ${_project}_DEFINITIONS ${_option_name})
	if(DEFINED ${_prefix}_INCLUDE_DIRS)
          list(APPEND ${_project}_EXTERNAL_LIB_INCLUDE_DIR ${${_prefix}_INCLUDE_DIRS})
          set(${_pkg_name}_INCLUDE_DIR ${${_prefix}_INCLUDE_DIRS})
	elseif(DEFINED ${_prefix}_INCLUDE_DIR)
          list(APPEND ${_project}_EXTERNAL_LIB_INCLUDE_DIR ${${_prefix}_INCLUDE_DIR})
          set(${_pkg_name}_INCLUDE_DIR ${${_prefix}_INCLUDE_DIR})
	elseif(DEFINED ${_prefix}_INCLUDE_PATH)
          list(APPEND ${_project}_EXTERNAL_LIB_INCLUDE_DIR ${${_prefix}_INCLUDE_PATH})
          set(${_pkg_name}_INCLUDE_DIR ${${_prefix}_INCLUDE_PATH})
	endif()
	list(APPEND ${_project}_EXTERNAL_LIBRARIES ${${_prefix}_LIBRARIES})
	set(${_pkg_name}_LIBRARIES ${${_prefix}_LIBRARIES})
	set(${_pkg_name} ON)
        string(TOUPPER ${PACKAGE} _u_package)
	list(APPEND ${_project}_OPTION_LIST ${_u_package})
	cmake_debug_message(PackagesSystem "add_optional_package: Package ${PACKAGE} found! (PREFIX: ${_prefix})")
	cmake_debug_message(PackagesSystem "add_optional_package: Package ${PACKAGE} includes : ${${_pkg_name}_INCLUDE_DIR}")
	cmake_debug_message(PackagesSystem "add_optional_package: Package ${PACKAGE} libraries: ${${_pkg_name}_LIBRARIES}")
	cmake_debug_message(PackagesSystem "add_optional_package: option list: ${${_project}_OPTION_LIST}")
      else(${_prefix}_FOUND)
	cmake_debug_message(PackagesSystem "add_optional_package: Package ${PACKAGE} not found! (PREFIX: ${_prefix})")
	set(${_pkg_name} OFF)
      endif()
    endforeach()
  endif(${_option_name})
endmacro()

#===============================================================================
# macro to add meta packages
macro(add_meta_package PKG DESC DEFAULT)
  cmake_debug_message(PackagesSystem "add_meta_package: register meta option ${PKG} ${DESC} ${DEFAULT}")
  package_pkg_name (${PKG} _pkg_name)
  package_desc_name(${PKG} _desc_name)

  set(${_desc_name} ${DESC})
  option(${_pkg_name} ${DESC} ${DEFAULT})

  foreach(_dep ${ARGN})
    package_opt_name (${_dep} _dep_name)
    mark_as_advanced(${_dep_name})
    add_external_package_dependencies(${PKG} ${_dep})
  endforeach()
endmacro()

#===============================================================================
macro(_add_package_dependencies PKG DEP _dep_name)
  package_pkg_name (${PKG} _opt_name)
  package_desc_name(${DEP} _var_dep_desc)

  cmake_debug_message(PackagesSystem "add_package_dependecies: add dependence between ${_opt_name} and ${_dep_name}")
  set(_dep_desc ${_var_dep_desc})

  cmake_debug_message(PackagesSystem "add_package_dependecies: ON dependencies of ${_dep_name} are: ${${_dep_name}_DEPS}")
  cmake_debug_message(PackagesSystem "add_package_dependecies: saved value for ${_dep_name} is: ${${_dep_name}_OLD}")
  if(${_opt_name})
    if("${${_dep_name}_DEPS}" STREQUAL "")
      cmake_debug_message(PackagesSystem "add_package_dependecies: Save dep state ${_dep_name}:${${_dep_name}}")
      set(${_dep_name}_OLD ${${_dep_name}} CACHE INTERNAL "${_dep_desc}" FORCE)
    endif()

    cmake_debug_message(PackagesSystem "add_package_dependecies: force value to ON ${_dep_name}")
    set(${_dep_name} ON CACHE BOOL "${_dep_desc}" FORCE)

    list(FIND ${_dep_name}_DEPS ${_opt_name} pos)
    if(pos EQUAL -1)
      list(APPEND ${_dep_name}_DEPS ${_opt_name})
      set(${_dep_name}_DEPS ${${_dep_name}_DEPS} CACHE INTERNAL "Dependencies ON with package ${_dep_name}" FORCE)
    endif()
  else()
    list(LENGTH ${_dep_name}_DEPS len)
    list(FIND ${_dep_name}_DEPS ${_opt_name} pos)
    if((len EQUAL 1) AND (NOT pos EQUAL -1))
      cmake_debug_message(PackagesSystem "add_package_dependecies: Restore state ${_dep_name}:${${_dep_name}} (${pos})")
      set(${_dep_name} ${${_dep_name}_OLD} CACHE BOOL "${_dep_desc}" FORCE)
      unset(${_dep_name}_OLD CACHE)
    endif()

    if(NOT pos EQUAL -1)
      list(REMOVE_AT ${_dep_name}_DEPS ${pos})
      set(${_dep_name}_DEPS ${${_dep_name}_DEPS} CACHE INTERNAL "Nb dependencies with package ${_dep_name}" FORCE)
    endif()
  endif()
endmacro()

#===============================================================================
macro(add_internal_package_dependencies PKG DEP)
  package_pkg_name (${DEP} _dep_name)
  _add_package_dependencies(${PKG} ${DEP} ${_dep_name})
endmacro()

#===============================================================================
macro(add_external_package_dependencies PKG DEP)
  package_opt_name (${DEP} _dep_name)
  _add_package_dependencies(${PKG} ${DEP} ${_dep_name})
endmacro()
