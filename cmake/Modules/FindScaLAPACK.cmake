#===============================================================================
# @file   FindMumps.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Dec 13 2010
# @date last modification: Tue Sep 09 2014
#
# @brief  The find_package file for the Mumps solver
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

find_library(BLACS_LIBRARY_C NAME blacsC
  HINTS ${SCALAPACK_DIR} PATH_SUFFIXES lib)

find_library(BLACS_LIBRARY_F77 NAME blacsF77
  HINTS ${SCALAPACK_DIR} PATH_SUFFIXES lib)

find_library(BLACS_LIBRARY NAME blacs
  HINTS ${SCALAPACK_DIR} PATH_SUFFIXES lib)

find_library(SCALAPACK_LIBRARY NAME scalapack
  HINTS ${SCALAPACK_DIR} PATH_SUFFIXES lib)

mark_as_advanced(BLACS_LIBRARY_C)
mark_as_advanced(BLACS_LIBRARY_F77)
mark_as_advanced(BLACS_LIBRARY)
mark_as_advanced(SCALAPACK_LIBRARY)
if(SCALAPACK_LIBRARY)
  list(APPEND SCALAPACK_LIBRARIES_ALL ${SCALAPACK_LIBRARY})
endif()
if(BLACS_LIBRARY_F77)
  list(APPEND SCALAPACK_LIBRARIES_ALL ${BLACS_LIBRARY_F77})
endif()
if(BLACS_LIBRARY)
  list(APPEND SCALAPACK_LIBRARIES_ALL ${BLACS_LIBRARY})
endif()
if(BLACS_LIBRARY_C)
  list(APPEND SCALAPACK_LIBRARIES_ALL ${BLACS_LIBRARY_C})
endif()
if(BLACS_LIBRARY_F77)
  list(APPEND SCALAPACK_LIBRARIES_ALL ${BLACS_LIBRARY_F77})
endif()

set(SCALAPACK_LIBRARIES ${SCALAPACK_LIBRARIES_ALL} CACHE INTERNAL "Libraries for MUMPS" FORCE)

#===============================================================================
if(NOT MUMPS_FOUND)
  set(SCALAPACK_DIR "" CACHE PATH "Prefix of MUMPS library.")
  mark_as_advanced(SCALAPACK_DIR)
endif()

#===============================================================================
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ScaLAPACK DEFAULT_MSG
  SCALAPACK_LIBRARIES SCALAPACK_INCLUDE_DIR)