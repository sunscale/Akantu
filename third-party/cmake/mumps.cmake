#===============================================================================
# @file   mumps.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Nov 21 2011
# @date last modification: Mon Sep 15 2014
#
# @brief  compilation of the third-party MUMPS
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

enable_language(Fortran)

set(MUMPS_VERSION ${AKANTU_USE_MUMPS_VERSION})
set(MUMPS_ARCHIVE "${PROJECT_SOURCE_DIR}/third-party/MUMPS_${MUMPS_VERSION}.tar.gz")
set(MUMPS_ARCHIVE_HASH_4.9.2  "MD5=d0b8f139a4acf29b76dbae69ade8ac54")
set(MUMPS_ARCHIVE_HASH_4.10.0 "MD5=959e9981b606cd574f713b8422ef0d9f")
set(MUMPS_ARCHIVE_HASH_5.0.0  "MD5=3c6aeab847e9d775ca160194a9db2b75")

if(NOT EXISTS ${MUMPS_ARCHIVE})
  message(FATAL_ERROR "To be able to compile MUMPS please download it from "
    "http://mumps.enseeiht.fr/ or http://graal.ens-lyon.fr/MUMPS and place it "
    "in the directory: ${PROJECT_SOURCE_DIR}/third-party and in cmake set the "
    "variable AKANTU_MUMPS_VERSION to the corresponding version "
    "\n"
    "Supported version for automated compilation in Akantu are 4.9.2, 4.10.0 "
    "and 5.0.0")
endif()

package_get_option_name(MPI _mpi_option)
if(${_mpi_option})
  unset(MUMPS_PREFIX)
else()
  set(MUMPS_PREFIX _seq)
endif()

package_use_system(Scotch _scotch_use_system)
if(NOT _scotch_use_system)
  list(APPEND MUMPS_DEPENDS Scotch)
endif()

package_get_option_name(ScaLAPACK _scalapack_option)
package_use_system(ScaLAPACK _scalapack_use_system)
if(NOT _scalapack_use_system AND ${_scalapack_option})
  list(APPEND MUMPS_DEPENDS ScaLAPACK)
endif()

include(AkantuMacros)

package_get_libraries(Scotch _scotch_libraries)
string(REPLACE ";" " " MUMPS_SCOTCH_LIBRARIES
  "${_scotch_libraries};${SCOTCH_LIBRARY_ESMUMPS}")

package_get_libraries(BLAS _blas_libraries)
foreach(_blas_lib ${_blas_libraries})
  if("${_blas_lib}" MATCHES ".*\\.framework")
    get_filename_component(_blas_framework "${_blas_lib}" NAME_WE)
    set(MUMPS_BLAS_LIBRARIES "${MUMPS_BLAS_LIBRARIES} -framework ${_blas_framework}")
  else()
    set(MUMPS_BLAS_LIBRARIES "${MUMPS_BLAS_LIBRARIES} ${_blas_lib}")
  endif()
endforeach()

if("${MUMPS_TYPE}" STREQUAL "seq")
  set_third_party_shared_libirary_name(MUMPS_LIBRARY_MPI mpiseq${MUMPS_PREFIX})
  mark_as_advanced(MUMPS_LIBRARY_MPI)
else()
  set(MUMPS_LIBRARY_MPI "")
endif()

if(CMAKE_C_COMPILER_ID STREQUAL "Intel")
  set(MUMPS_EXTRA_Fortran_FLAGS "-nofor_main")
else()
  set(MUMPS_EXTRA_Fortran_FLAGS "")
endif()

if(CMAKE_VERSION VERSION_GREATER 3.1)
  set(_extra_options
    DOWNLOAD_NO_PROGRESS 1
    EXCLUDE_FROM_ALL 1
    )
endif()

if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
  set(AKANTU_MUMPS_CDEFS "-DAdd_ -DWITHOUT_PTHREAD")
else()
  set(AKANTU_MUMPS_CDEFS "-DAdd_")
  set(AKANTU_MUMPS_PTHREAD "-lpthread")
endif()

configure_file(${PROJECT_SOURCE_DIR}/third-party/MUMPS_${MUMPS_VERSION}_make.inc.cmake
  ${PROJECT_BINARY_DIR}/third-party/MUMPS_make.inc @ONLY)

ExternalProject_Add(MUMPS
  DEPENDS ${MUMPS_DEPENDS}
  PREFIX ${PROJECT_BINARY_DIR}/third-party
  URL ${MUMPS_ARCHIVE}
  URL_HASH ${MUMPS_ARCHIVE_HASH_${MUMPS_VERSION}}
  ${_extra_options}
  BUILD_IN_SOURCE 1
  PATCH_COMMAND ${PATCH_COMMAND} -p2 < ${PROJECT_SOURCE_DIR}/third-party/MUMPS_${MUMPS_VERSION}.patch
  CONFIGURE_COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_BINARY_DIR}/third-party/MUMPS_make.inc Makefile.inc
  BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} d
  INSTALL_COMMAND "${CMAKE_MAKE_PROGRAM}" prefix=<INSTALL_DIR> install
  LOG_DOWNLOAD 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
  )

set_third_party_shared_libirary_name(MUMPS_LIBRARY_DMUMPS dmumps${MUMPS_PREFIX})
set_third_party_shared_libirary_name(MUMPS_LIBRARY_COMMON mumps_common${MUMPS_PREFIX})
set_third_party_shared_libirary_name(MUMPS_LIBRARY_PORD pord${MUMPS_PREFIX})

mark_as_advanced(
  MUMPS_LIBRARY_COMMON
  MUMPS_LIBRARY_DMUMPS
  MUMPS_LIBRARY_PORD
  MUMPS_LIBRARY_MPI
  MUMPS_INCLUDE_DIR
  )

set(MUMPS_LIBRARIES_ALL
  ${MPI_Fortran_LIBRARIES}
  ${MUMPS_LIBRARY_COMMON}
  ${MUMPS_LIBRARY_DMUMPS}
  ${MUMPS_LIBRARY_PORD}
  ${MUMPS_LIBRARY_MPI})

set(MUMPS_INCLUDE_DIR ${PROJECT_BINARY_DIR}/third-party/include CACHE PATH "" FORCE)
set(MUMPS_LIBRARIES ${MUMPS_LIBRARIES_ALL} CACHE INTERNAL "Libraries for MUMPS" FORCE)

package_add_extra_dependency(Mumps MUMPS)
