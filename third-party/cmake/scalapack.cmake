#===============================================================================
# @file   scalapack.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Mar 30 2015
# @date last modification: Wed Jun 10 2015
#
# @brief  package description for mumps support
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

enable_language(Fortran)

if(CMAKE_VERSION VERSION_GREATER 3.1)
  set(_extra_options 
    UPDATE_DISCONNECTED 1
    DOWNLOAD_NO_PROGRESS 1
    EXCLUDE_FROM_ALL 1
    )
endif()

set(SCALAPACK_DIR ${PROJECT_BINARY_DIR}/third-party)
ExternalProject_Add(ScaLAPACK
  PREFIX   ${SCALAPACK_DIR}
  URL      ${SCALAPACK_ARCHIVE}
  URL_HASH ${SCALAPACK_ARCHIVE_HASH}
  ${_extra_options}
  PATCH_COMMAND patch -p1 < ${PROJECT_SOURCE_DIR}/third-party/scalapack_${SCALAPACK_VERSION}.patch
  CMAKE_COMMAND BLA_VENDOR=$ENV{BLA_VENDOR} MPICH_F90=${CMAKE_Fortran_COMPILER} ${CMAKE_COMMAND}
  CMAKE_ARGS <SOURCE_DIR>/ScaLAPACK
  CMAKE_CACHE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
  -DMPIEXEC:PATH=${MPIEXEC}
  -DCMAKE_C_COMPILER:PATH=${CMAKE_C_COMPILER}
  -DMPI_C_COMPILER:PATH=${MPI_C_COMPILER}
  -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
  -DCMAKE_Fortran_COMPILER:PATH=${CMAKE_Fortran_COMPILER}
  -DCMAKE_Fortran_FLAGS:STRING=${CMAKE_Fortran_FLAGS}
  -DMPI_Fortran_COMPILER:PATH=${MPI_Fortran_COMPILER}
  -DBUILD_SHARED_LIBS:BOOL=ON
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  BUILD_COMMAND BLA_VENDOR=$ENV{BLA_VENDOR} MPICH_F90=${CMAKE_Fortran_COMPILER} ${CMAKE_MAKE_PROGRAM}
  LOG_DOWNLOAD 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
  )

set_third_party_shared_libirary_name(SCALAPACK_LIBRARIES scalapack)
set(SCALAPACK_INCLUDE_DIR ${SCALAPACK_DIR}/include CACHE PATH "" FORCE)
mark_as_advanced(SCALAPACK_LIBRARIES SCALAPACK_INCLUDE_DIR)
