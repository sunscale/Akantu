#===============================================================================
# @file   CMakePackagesSystemGlobalFunctions.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Sat Jul 18 2015
# @date last modification: Mon Jan 18 2016
#
# @brief  Set of macros used by the package system to set internal variables
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

# ==============================================================================
# Package system meta functions
# ==============================================================================
function(package_set_project_variable variable)
  string(TOUPPER ${PROJECT_NAME} _u_project)
  set(${_u_project}_${variable} "${ARGN}" CACHE INTERNAL "" FORCE)
endfunction()

function(package_get_project_variable variable value_out)
  string(TOUPPER ${PROJECT_NAME} _u_project)
  set(${value_out} ${${_u_project}_${variable}} PARENT_SCOPE)
endfunction()

function(package_add_to_project_variable variable)
  package_get_project_variable(${variable} _tmp_list)
  list(APPEND _tmp_list ${ARGN})
  if(_tmp_list)
    list(REMOVE_DUPLICATES _tmp_list)
  endif()
  package_set_project_variable(${variable} ${_tmp_list})
endfunction()

# ==============================================================================
function(_package_set_variable variable pkg_name)
  set(${pkg_name}_${variable} ${ARGN} CACHE INTERNAL "" FORCE)
endfunction()

function(_package_get_variable variable pkg_name value)
  #unset(${value} PARENT_SCOPE)
  if(DEFINED ${pkg_name}_${variable})
    set(${value} ${${pkg_name}_${variable}} PARENT_SCOPE)
  elseif(DEFINED ARGN)
    set(${value} ${ARGN} PARENT_SCOPE)
  else()
    set(${value} PARENT_SCOPE)
  endif()
endfunction()

# ==============================================================================
function(_package_variable_unset variable pkg_name)
  unset(${pkg_name}_${variable} CACHE)
endfunction()

# ==============================================================================
function(_package_add_to_variable variable pkg_name)
  _package_get_variable(${variable} ${pkg_name} _tmp_list)
  list(APPEND _tmp_list ${ARGN})
  if(_tmp_list)
    list(REMOVE_DUPLICATES _tmp_list)
  endif()
  _package_set_variable(${variable} ${pkg_name} ${_tmp_list})
endfunction()

function(_package_remove_from_variable variable pkg_name value)
  _package_get_variable(${variable} ${pkg_name} _tmp_list)
  list(LENGTH _tmp_list _length)
  if(_length GREATER 0)
    list(REMOVE_ITEM _tmp_list ${value})
    _package_set_variable(${variable} ${pkg_name} ${_tmp_list})
  endif()
endfunction()

# ==============================================================================
function(_package_get_variable_for_packages variable values)
  set(_list_values)
  foreach(_pkg_name ${ARGN})
    _package_get_variable(${variable} ${_pkg_name} _value)
    list(APPEND _list_values ${_value})
  endforeach()

  if (_list_values)
    list(REMOVE_DUPLICATES _list_values)
  endif()

  set(${values} ${_list_values} PARENT_SCOPE)
endfunction()

# ==============================================================================
function(_package_get_variable_for_activated variable values)
  package_get_all_activated_packages(_activated_list)
  _package_get_variable_for_packages(${variable} _list_values ${_activated_list})
  set(${values} ${_list_values} PARENT_SCOPE)
endfunction()

# ==============================================================================
function(_package_get_variable_for_external_dependencies variable type values)
  set(_list_packages)
  package_get_all_activated_packages(_activated_list)
  foreach(_pkg_name ${_activated_list})
    _package_get_nature(${_pkg_name} _nature)
    if(_nature MATCHES "^external.*")
      list(APPEND _list_packages ${_pkg_name})
    endif()
  endforeach()

  if (_list_packages)
    list(REMOVE_DUPLICATES _list_packages)
    _package_get_variable_for_packages(${variable} _list_values_deps ${_list_packages})
    set(${values} ${_list_values_deps} PARENT_SCOPE)
  endif()
endfunction()
# ==============================================================================
