#===============================================================================
# @file   CMakePackagesSystem.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @brief  Addition to the PackageSystem specific for Akantu
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

#===============================================================================
# Material specific
#===============================================================================
#-------------------------------------------------------------------------------
function(package_declare_material_infos pkg)
  cmake_parse_arguments(_opt_pkg
    ""
    ""
    "LIST;INCLUDE"
    ${ARGN})

  package_set_variable(MATERIAL_LIST ${pkg} ${_opt_pkg_LIST})
  package_set_variable(MATERIAL_INCLUDE ${pkg} ${_opt_pkg_INCLUDE})
endfunction()

#-------------------------------------------------------------------------------
function(package_get_all_material_includes includes)
  _package_get_variable_for_activated(MATERIAL_INCLUDE _includes)

  foreach(_mat_inc ${_includes})
    if(DEFINED _mat_includes)
      set(_mat_includes "${_mat_includes}\n#include \"${_mat_inc}\"")
    else()
      set(_mat_includes "#include \"${_mat_inc}\"")
    endif()
  endforeach()

  set(${includes} ${_mat_includes} PARENT_SCOPE)
endfunction()

#-------------------------------------------------------------------------------
function(package_get_all_material_lists lists)
  _package_get_variable_for_activated(MATERIAL_LIST _lists)

  foreach(_mat_list ${_lists})
    if(DEFINED _mat_lists)
      set(_mat_lists "${_mat_lists}\n  ${_mat_list}\t\t\t\\")
    else()
      set(_mat_lists "  ${_mat_list}\t\t\t\\")
    endif()
  endforeach()

  set(${lists} ${_mat_lists} PARENT_SCOPE)
endfunction()

#-------------------------------------------------------------------------------

