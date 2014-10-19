#===============================================================================
# @file   85_mumps.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Nov 21 2011
# @date last modification: Mon Sep 15 2014
#
# @brief  package description for mumps support
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

set(AKANTU_MUMPS_FILES
  solver/solver_mumps.cc
  solver/solver_mumps.hh
  )

option(AKANTU_USE_THIRD_PARTY_MUMPS "Use the third-party Mumps instead of the one from the system" OFF)
option(AKANTU_USE_MUMPS "Add Mumps support in akantu")


mark_as_advanced(AKANTU_USE_THIRD_PARTY_MUMPS)
if(AKANTU_USE_THIRD_PARTY_MUMPS AND AKANTU_USE_MUMPS)
  set(MUMPS_DEPENDS)
  enable_language(Fortran)

  include(AkantuMacros)
  include(ExternalProject)
  if(AKANTU_USE_MPI)
    set(SCALAPACK_VERSION "2.0.2")
    set(SCALAPACK_ARCHIVE "http://www.netlib.org/scalapack/scalapack-${SCALAPACK_VERSION}.tgz")
    set(SCALAPACK_ARCHIVE_HASH "MD5=2f75e600a2ba155ed9ce974a1c4b536f")

    ExternalProject_Add(ScaLAPACK
      PREFIX ${PROJECT_BINARY_DIR}/third-party
      URL      ${SCALAPACK_ARCHIVE}
      URL_HASH ${SCALAPACK_ARCHIVE_HASH}
      PATCH_COMMAND patch -p1 < ${PROJECT_SOURCE_DIR}/third-party/scalapack_${SCALAPACK_VERSION}.patch
      CMAKE_ARGS <SOURCE_DIR>/ScaLAPACK
      CMAKE_CACHE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> -DCMAKE_C_FLAGS:STRING=-fPIC -DCMAKE_Fortran_FLAGS:STRING=-fPIC -DBUILD_SHARED_LIBS:BOOL=ON
      )

    set_third_party_shared_libirary_name(SCALAPACK_LIBRARIES scalapack)
    mark_as_advanced(SCALAPACK_LIBRARIES)

    install(FILES ${SCALAPACK_LIBRARIES}
      DESTINATION lib
      COMPONENT lib)

    list(APPEND MUMPS_DEPENDS ScaLAPACK)
    set(MUMPS_TYPE par)
    unset(MUMPS_PREFIX)
  else()
    set(MUMPS_TYPE seq)
    set(MUMPS_PREFIX _seq)
  endif()

  if(AKANTU_USE_THIRD_PARTY_SCOTCH)
    include(${PROJECT_SOURCE_DIR}/packages/90_scotch.cmake)
    list(APPEND MUMPS_DEPENDS Scotch)
  else()
    find_package(Scotch REQUIRED)
  endif()
  string(REPLACE ";" " " MUMPS_SCOTCH_LIBRARIES "${SCOTCH_LIBRARIES};${SCOTCH_LIBRARY_ESMUMPS}")

  find_package(BLAS REQUIRED)
  foreach(_blas_lib ${BLAS_LIBRARIES})
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

  set(MUMPS_VERSION "4.10.0")
  set(MUMPS_ARCHIVE "${PROJECT_SOURCE_DIR}/third-party/MUMPS_${MUMPS_VERSION}.tar.gz")
  set(MUMPS_ARCHIVE_HASH "MD5=959e9981b606cd574f713b8422ef0d9f")
  if(NOT EXISTS ${MUMPS_ARCHIVE})
    set(MUMPS_VERSION "4.9.2")
    set(MUMPS_ARCHIVE "${PROJECT_SOURCE_DIR}/third-party/MUMPS_${MUMPS_VERSION}.tar.gz")
    set(MUMPS_ARCHIVE_HASH "MD5=d0b8f139a4acf29b76dbae69ade8ac54")
    if(NOT EXISTS ${MUMPS_ARCHIVE})
      MESSAGE(ERROR "To be able to compile MUMPS please download it from
http://mumps.enseeiht.fr/ or http://graal.ens-lyon.fr/MUMPS
and place it in the directory: ${PROJECT_SOURCE_DIR}/third-party

Supported version for automated compilation in Akantu are 4.9.2 and 4.10.0")
    endif()
  endif()

  configure_file(${PROJECT_SOURCE_DIR}/third-party/MUMPS_${MUMPS_VERSION}_make.inc.cmake
    ${PROJECT_BINARY_DIR}/third-party/MUMPSmake.inc)

  ExternalProject_Add(MUMPS
    DEPENDS ${MUMPS_DEPENDS}
    PREFIX ${PROJECT_BINARY_DIR}/third-party
    URL ${MUMPS_ARCHIVE}
    URL_HASH ${MUMPS_ARCHIVE_HASH}
    BUILD_IN_SOURCE 1
    PATCH_COMMAND patch -p2 < ${PROJECT_SOURCE_DIR}/third-party/MUMPS_${MUMPS_VERSION}.patch
    CONFIGURE_COMMAND cmake -E copy ${PROJECT_BINARY_DIR}/third-party/MUMPSmake.inc Makefile.inc
    BUILD_COMMAND make d
    INSTALL_COMMAND prefix=<INSTALL_DIR> make install
    )

  set_third_party_shared_libirary_name(MUMPS_LIBRARY_DMUMPS dmumps${MUMPS_PREFIX})
  set_third_party_shared_libirary_name(MUMPS_LIBRARY_COMMON mumps_common${MUMPS_PREFIX})
  set_third_party_shared_libirary_name(MUMPS_LIBRARY_PORD pord${MUMPS_PREFIX})
  set(MUMPS_INCLUDE_DIR ${PROJECT_BINARY_DIR}/third-party/include CACHE PATH "" FORCE)
  set(MUMPS_LIBRARIES_ALL)

  mark_as_advanced(MUMPS_LIBRARY_COMMON MUMPS_LIBRARY_DMUMPS MUMPS_LIBRARY_PORD MUMPS_INCLUDE_DIR)

  list(APPEND MUMPS_LIBRARIES_ALL ${MPI_Fortran_LIBRARIES} ${MUMPS_LIBRARY_COMMON} ${MUMPS_LIBRARY_DMUMPS} ${MUMPS_LIBRARY_PORD} ${MUMPS_LIBRARY_MPI})
  set(MUMPS_LIBRARIES ${MUMPS_LIBRARIES_ALL} CACHE INTERNAL "Libraries for MUMPS" FORCE)

  install(FILES ${MUMPS_LIBRARIES}
    DESTINATION lib
    COMPONENT lib)
  install(DIRECTORY ${PROJECT_BINARY_DIR}/third-party/include/ DESTINATION include/mumps
    COMPONENT dev
    FILES_MATCHING PATTERN "*mumps*.h")

  list(APPEND AKANTU_EXTERNAL_LIBRARIES ${MUMPS_LIBRARIES})
  list(APPEND AKANTU_EXTERNAL_LIB_INCLUDE_DIR ${MUMPS_INCLUDE_DIR})

  set(AKANTU_MUMPS_INCLUDE_DIR ${MUMPS_INCLUDE_DIR})
  set(AKANTU_MUMPS_LIBRARIES ${MUMPS_LIBRARIES})

  list(APPEND AKANTU_OPTION_LIST MUMPS)

  set(MUMPS_FOUND TRUE CACHE INTERNAL "" FORCE)
  set(AKANTU_MUMPS ON)

  list(APPEND AKANTU_EXTRA_TARGET_DEPENDENCIES MUMPS)
else()
  if(AKANTU_USE_MPI)
    set(MUMPS_TYPE par)
    set(AKANTU_MUMPS_DEB_DEPEND
      libmumps-dev
      )
  else()
    set(MUMPS_TYPE seq)
    set(AKANTU_MUMPS_DEB_DEPEND
      libmumps-seq-dev
      )
  endif()

  add_optional_external_package(Mumps "Add Mumps support in akantu" OFF)

endif()

set(AKANTU_MUMPS_TESTS
  test_sparse_matrix_profile
  test_sparse_matrix_assemble
  test_solver_mumps
  test_sparse_matrix_product
  )

set(AKANTU_MUMPS_DOCUMENTATION "
This package enables the \\href{http://mumps.enseeiht.fr/}{MUMPS} parallel direct solver for sparce matrices.
This is necessary to solve static or implicit problems.

Under Ubuntu (14.04 LTS) the installation can be performed using the commands:

\\begin{command}
  > sudo apt-get install libmumps-seq-dev # for sequential
  > sudo apt-get install libmumps-dev     # for parallel
\\end{command}

Under Mac OS X the installation requires the following steps:
\\begin{command}
  > sudo port install mumps
\\end{command}

If you activate the advanced option AKANTU\\_USE\\_THIRD\\_PARTY\\_MUMPS the make system of akantu can automatically compile MUMPS. For this you will have to download MUMPS from \\url{http://mumps.enseeiht.fr/} or \\url{http://graal.ens-lyon.fr/MUMPS} and place it in \\shellcode{<akantu source>/third-party}
")
