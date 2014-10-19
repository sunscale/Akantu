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

macro(include_package_if_needed PACKAGE)
  string(TOUPPER ${PACKAGE} _package)
  if(AKANTU_USE_${_package})
    if(${PACKAGE} MATCHES BLAS)
      enable_language(Fortran)
    endif()
    find_package(${PACKAGE} REQUIRED)
    if(${_package}_FOUND)
      if(DEFINED ${_package}_INCLUDE_DIR)
        list(APPEND AKANTU_EXTRA_INCLUDE_DIR ${${_package}_INCLUDE_DIR})
      else()
        list(APPEND AKANTU_EXTRA_INCLUDE_DIR ${${_package}_INCLUDE_PATH})
      endif()
      list(APPEND AKANTU_EXTRA_LIBRARIES ${${_package}_LIBRARIES})
    endif()
  endif()
endmacro()

macro(find_akantu_dependencies)
  include_package_if_needed(MPI)
  include_package_if_needed(Mumps)
  include_package_if_needed(Scotch)
  include_package_if_needed(BLAS)
endmacro()