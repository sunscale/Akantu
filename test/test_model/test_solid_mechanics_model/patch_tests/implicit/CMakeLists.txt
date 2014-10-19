#===============================================================================
# @file   CMakeLists.txt
#
# @author David Simon Kammer <david.kammer@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Thu Feb 17 2011
# @date last modification: Thu Feb 21 2013
#
# @brief  configuration for the implicit patch test
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
# @section DESCRIPTION
#
#===============================================================================

macro(register_patch_test name type plane_strain)
  if("${plane_strain}" STREQUAL "false")
    set(_mat ../data/material_check_stress_plane_stress.dat)
  else()
    set(_mat ../data/material_check_stress_plane_strain.dat)
  endif()
  register_test(patch_test_implicit_${name}
    SOURCES patch_test_implicit.cc
    FILES_TO_COPY ../data/_${type}.msh ${_mat}
    COMPILE_OPTIONS "TYPE=_${type};PLANE_STRAIN=${plane_strain}"
    PACKAGE implicit
    )
endmacro()

foreach(_type ${LIST_TYPES})
  register_patch_test(${_type} ${_type} true)
endforeach()

foreach(_type ${LIST_TYPES_2D})
  register_patch_test(plane_strain_${_type} ${_type} true)
  register_patch_test(plane_stress_${_type} ${_type} false)
endforeach()
