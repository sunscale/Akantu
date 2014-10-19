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

#===============================================================================
if(NOT MUMPS_TYPE)
  set(MUMPS_TYPE par)
endif()

if("${MUMPS_TYPE}" STREQUAL "seq")
  set(MUMPS_PREFIX _seq)
else()
  unset(MUMPS_PREFIX)
endif()

find_library(MUMPS_LIBRARY_DMUMPS NAMES dmumps${MUMPS_PREFIX}
   HINTS ${MUMPS_DIR} /usr
   PATH_SUFFIXES lib
   )

find_library(MUMPS_LIBRARY_COMMON NAMES mumps_common${MUMPS_PREFIX}
   HINTS ${MUMPS_DIR}
   PATH_SUFFIXES lib
   )

find_library(MUMPS_LIBRARY_PORD NAMES pord${MUMPS_PREFIX}
   HINTS ${MUMPS_DIR}
   PATH_SUFFIXES lib
   )


find_path(MUMPS_INCLUDE_DIR dmumps_c.h
  HINTS ${MUMPS_DIR}
  PATH_SUFFIXES include
  )

mark_as_advanced(MUMPS_LIBRARY_COMMON)
mark_as_advanced(MUMPS_LIBRARY_DMUMPS)
mark_as_advanced(MUMPS_LIBRARY_PORD)
mark_as_advanced(MUMPS_INCLUDE_DIR)
set(MUMPS_LIBRARIES_ALL ${MUMPS_LIBRARY_DMUMPS} ${MUMPS_LIBRARY_COMMON} ${MUMPS_LIBRARY_PORD})

if("${MUMPS_TYPE}" STREQUAL "par")
  find_library(BLACS_LIBRARY_C NAME blacsC
    HINTS ${MUMPS_DIR} PATH_SUFFIXES lib)
  find_library(BLACS_LIBRARY_F77 NAME blacsF77
    HINTS ${MUMPS_DIR} PATH_SUFFIXES lib)
  find_library(BLACS_LIBRARY NAME blacs
    HINTS ${MUMPS_DIR} PATH_SUFFIXES lib)
  find_library(SCALAPACK_LIBRARIES NAME scalapack
    HINTS ${MUMPS_DIR} PATH_SUFFIXES lib)

  mark_as_advanced(BLACS_LIBRARY_C)
  mark_as_advanced(BLACS_LIBRARY_F77)
  mark_as_advanced(BLACS_LIBRARY)
  mark_as_advanced(SCALAPACK_LIBRARY)
  mark_as_advanced(SCALAPACK_LIBRARIES)
  if(SCALAPACK_LIBRARY)
    set(BLACS_LIBRARIES_ALL ${BLACS_LIBRARIES_ALL} ${SCALAPACK_LIBRARY})
  endif()
  if(BLACS_LIBRARY_F77)
    set(BLACS_LIBRARIES_ALL ${BLACS_LIBRARIES_ALL} ${BLACS_LIBRARY_F77})
  endif()
  if(BLACS_LIBRARY)
    set(BLACS_LIBRARIES_ALL ${BLACS_LIBRARIES_ALL} ${BLACS_LIBRARY})
  endif()
  if(BLACS_LIBRARY_C)
    set(BLACS_LIBRARIES_ALL ${BLACS_LIBRARIES_ALL} ${BLACS_LIBRARY_C})
  endif()
  if(BLACS_LIBRARY_F77)
    set(BLACS_LIBRARIES_ALL ${BLACS_LIBRARIES_ALL} ${BLACS_LIBRARY_F77})
  endif()
endif()

set(MUMPS_LIBRARIES ${MUMPS_LIBRARIES_ALL} ${BLACS_LIBRARIES_ALL} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}  CACHE INTERNAL "Libraries for MUMPS" FORCE)


#===============================================================================
if(NOT MUMPS_FOUND)
  set(MUMPS_DIR "" CACHE PATH "Prefix of MUMPS library.")
  mark_as_advanced(MUMPS_DIR)
endif()
#===============================================================================
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Mumps DEFAULT_MSG
  MUMPS_LIBRARIES MUMPS_INCLUDE_DIR)