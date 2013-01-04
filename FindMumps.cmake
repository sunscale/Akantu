#===============================================================================
# @file   FindMumps.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date   Mon Dec 13 10:48:06 2010
#
# @brief  The find_package file for the Mumps solver
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
find_library(MUMPS_LIBRARY_DMUMPS NAMES dmumps_seq dmumps_ptscotch dmumps
   PATHS ${MUMPS_DIR} /usr
   PATH_SUFFIXES lib
   )

find_library(MUMPS_LIBRARY_COMMON NAMES mumps_seq mumps_common_ptscotch mumps_common
   PATHS ${MUMPS_DIR}
   PATH_SUFFIXES lib
   )

find_library(MUMPS_LIBRARY_PORD NAMES pord_seq pord_ptscotch pord
   PATHS ${MUMPS_DIR}
   PATH_SUFFIXES lib
   )


find_path(MUMPS_INCLUDE_DIR dmumps_c.h
  PATHS ${MUMPS_DIR}
  PATH_SUFFIXES include
  )


find_library(BLACS_LIBRARY_C NAME blacsC
   PATHS ${MUMPS_DIR} PATH_SUFFIXES lib)
find_library(BLACS_LIBRARY_F77 NAME blacsF77
   PATHS ${MUMPS_DIR} PATH_SUFFIXES lib)
find_library(BLACS_LIBRARY NAME blacs
   PATHS ${MUMPS_DIR} PATH_SUFFIXES lib)
find_library(SCALAPACK_LIBRARY NAME scalapack
   PATHS ${MUMPS_DIR} PATH_SUFFIXES lib)

#enable_language(Fortran)
#find_package(BLAS REQUIRED)
#find_package(LAPACK REQUIRED)
##===============================================================================
mark_as_advanced(MUMPS_LIBRARY_COMMON)
mark_as_advanced(MUMPS_LIBRARY_DMUMPS)
mark_as_advanced(MUMPS_LIBRARY_PORD)
mark_as_advanced(MUMPS_INCLUDE_DIR)

mark_as_advanced(BLACS_LIBRARY_C)
mark_as_advanced(BLACS_LIBRARY_F77)
mark_as_advanced(BLACS_LIBRARY)
mark_as_advanced(SCALAPACK_LIBRARY)

set(MUMPS_LIBRARIES_ALL ${MUMPS_LIBRARY_DMUMPS} ${MUMPS_LIBRARY_COMMON} ${MUMPS_LIBRARY_PORD})

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

set(MUMPS_LIBRARIES ${MUMPS_LIBRARIES_ALL} ${BLACS_LIBRARIES_ALL} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}  CACHE INTERNAL "Libraries for MUMPS" FORCE)


#===============================================================================
if(NOT Mumps_FOUND)
  set(MUMPS_DIR "" CACHE PATH "Prefix of MUMPS library.")
endif(NOT Mumps_FOUND)
#===============================================================================
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Mumps DEFAULT_MSG
  MUMPS_LIBRARIES MUMPS_INCLUDE_DIR)