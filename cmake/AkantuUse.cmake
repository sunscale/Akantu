#===============================================================================
# @file   AkantuUse.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Thu Dec 01 2011
# @date last modification: Tue Nov 06 2012
#
# @brief  CMake file for the library
#
# @section LICENSE
#
# Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

function(package_is_activated pgk activated)
  string(TOUPPER ${pkg} _u_pkg)
  set(activated ${AKANTU_HAS_${_u_pkg}} PARENT_SCOPE)
endfunction()

function(package_get_include_dir pkg include_dir)
  string(TOUPPER ${pkg} _u_pkg)
  set(include_dir ${AKANTU_${_u_pkg}_INCLUDE_DIR} PARENT_SCOPE)
endfunction()

function(package_get_libraries pkg libs)
  string(TOUPPER ${pkg} _u_pkg)
  set(libs ${AKANTU_${_u_pkg}_LIBRARIES} PARENT_SCOPE)
endfunction()

function(package_get_compile_flags pkg flags)
  string(TOUPPER ${pkg} _u_pkg)
  set(flags ${AKANTU_${_u_pkg}_COMPILE_FLAGS} PARENT_SCOPE)
endfunction()
