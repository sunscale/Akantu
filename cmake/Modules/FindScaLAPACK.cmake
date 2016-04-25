#===============================================================================
# @file   FindScaLAPACK.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Tue Mar 31 2015
# @date last modification: Wed Jan 13 2016
#
# @brief  The find_package file for the Mumps solver
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

find_library(SCALAPACK_LIBRARY NAME scalapack
  PATHS "${SCALAPACK_DIR}" ENV SCALAPACK_DIR PATH_SUFFIXES lib)

mark_as_advanced(SCALAPACK_LIBRARY)

#===============================================================================
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ScaLAPACK DEFAULT_MSG
  SCALAPACK_LIBRARY)

if(SCALAPACK_FOUND AND NOT TARGET ScaLAPACK)
  set(SCALAPACK_LIBRARIES_ALL ${SCALAPACK_LIBRARY})

  include(CheckFortranFunctionExists)
  set(CMAKE_REQUIRED_LIBRARIES ${SCALAPACK_LIBRARIES_ALL})
  check_fortran_function_exists("blacs_gridinit" SCALAPACK_DOES_NOT_NEED_BLACS)
  set(CMAKE_REQUIRED_LIBRARIES)


  set(_blacs_dep)
  if(NOT SCALAPACK_DOES_NOT_NEED_BLACS)
    # Assuming scalapack was compiled as a static library
    set(SCALAPACK_LIBRARY_TYPE STATIC CACHE INTERNAL "" FORCE)

    find_library(BLACS_LIBRARY_C NAME blacsC
      PATHS "${SCALAPACK_DIR}" ENV SCALAPACK_DIR PATH_SUFFIXES lib)
    find_library(BLACS_LIBRARY_F77 NAME blacsF77
      PATHS "${SCALAPACK_DIR}" ENV SCALAPACK_DIR PATH_SUFFIXES lib)
    find_library(BLACS_LIBRARY NAME blacs
      PATHS "${SCALAPACK_DIR}" ENV SCALAPACK_DIR PATH_SUFFIXES lib)

    mark_as_advanced(
      BLACS_LIBRARY_C
      BLACS_LIBRARY_F77
      BLACS_LIBRARY
      )

    find_package_handle_standard_args(BLACS DEFAULT_MSG
      BLACS_LIBRARY BLACS_LIBRARY_C BLACS_LIBRARY_F77)

    add_library(blacs::common ${SCALAPACK_LIBRARY_TYPE} IMPORTED GLOBAL)
    add_library(blacs::F77 ${SCALAPACK_LIBRARY_TYPE} IMPORTED GLOBAL)
    add_library(blacs::C ${SCALAPACK_LIBRARY_TYPE} IMPORTED GLOBAL)

    set_target_properties(blacs::F77 PROPERTIES
      IMPORTED_LOCATION                 "${BLACS_LIBRARY_F77}"
      IMPORTED_LINK_INTERFACE_LANGUAGES "Fortran"
      INTERFACE_LINK_LIBRARIES  blacs::common
      )
    set_target_properties(blacs::C PROPERTIES
      IMPORTED_LOCATION                 "${BLACS_LIBRARY_C}"
      INTERFACE_INCLUDE_DIRECTORIES     "${SCALAPACK_INCLUDE_DIR}"
      IMPORTED_LINK_INTERFACE_LANGUAGES "C"
      INTERFACE_LINK_LIBRARIES          blacs::common
      )

    set_target_properties(blacs::common PROPERTIES
      IMPORTED_LOCATION                 "${BLACS_LIBRARY}"
      IMPORTED_LINK_INTERFACE_LANGUAGES "C Fortran"
      INTERFACE_LINK_LIBRARIES "blacs::C;blacs::F77"
      )

    find_package(LAPACK REQUIRED)
    find_package(BLAS REQUIRED)

    list(APPEND SCALAPACK_LIBRARIES_ALL
      ${BLACS_LIBRARY}
      ${BLACS_LIBRARY_C}
      ${BLACS_LIBRARY_F77}
      ${BLACS_LIBRARY}
      ${BLAS_LIBRARIES}
      ${LAPACK_LIBRARIES})

    set(_blacs_dep "blacs::common;${BLAS_LIBRARIES};${LAPACK_LIBRAIES}")
  else()
    set(SCALAPACK_LIBRARY_TYPE SHARED)
  endif()

  add_library(ScaLAPACK ${SCALAPACK_LIBRARY_TYPE} IMPORTED GLOBAL)
  set_target_properties(ScaLAPACK PROPERTIES
    IMPORTED_LOCATION                 "${SCALAPACK_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES     "${SCALAPACK_INCLUDE_DIR}"
    IMPORTED_LINK_INTERFACE_LANGUAGES "C Fortran"
    INTERFACE_LINK_LIBRARIES          "${_blacs_dep}")

  set(SCALAPACK_LIBRARIES ScaLAPACK CACHE INTERNAL "Libraries for ScaLAPACK" FORCE)
endif()
