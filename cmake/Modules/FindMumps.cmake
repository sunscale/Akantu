#===============================================================================
# @file   FindMumps.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Oct 24 2014
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

set(_MUMPS_COMPONENTS "sequential" "parallel" "double" "float" "complex_double" "complex_float")

if(NOT Mumps_FIND_COMPONENTS)
  set(Mumps_FIND_COMPONENTS "parallel" "double" "float" "complex_double" "complex_float")
endif()
#===============================================================================
enable_language(Fortran)

set(MUMPS_PRECISIONS)
foreach(_comp ${Mumps_FIND_COMPONENTS})
  if("${_comp}" STREQUAL "sequential")
    set(MUMPS_PLAT _seq) #default plat on debian based distribution
  endif()

  if("${_comp}" STREQUAL "float")
    list(APPEND MUMPS_PRECISIONS s)
  endif()
  if("${_comp}" STREQUAL "double")
    list(APPEND MUMPS_PRECISIONS d)
  endif()
  if("${_comp}" STREQUAL "complex_float")
    list(APPEND MUMPS_PRECISIONS c)
  endif()
  if("${_comp}" STREQUAL "complex_double")
    list(APPEND MUMPS_PRECISIONS z)
  endif()
endforeach()

if(NOT MUMPS_PRECISIONS)
  set(MUMPS_PRECISIONS s d c z)
endif()

list(GET MUMPS_PRECISIONS 0 _first_precision)

string(TOUPPER "${_first_precision}" _u_first_precision)

set(_mumps_test_dir "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}")
file(WRITE "${_mumps_test_dir}/mumps_test_code.c"
  "#include <${_first_precision}mumps_c.h>

int main() {
  ${_u_first_precision}MUMPS_STRUC_C mumps_data;
  ${_first_precision}mumps_c(&mumps_data);
  return 0;
}
")

find_path(MUMPS_INCLUDE_DIR ${_first_precision}mumps_c.h
  PATHS "${MUMPS_DIR}"
  ENV MUMPS_DIR
  PATH_SUFFIXES include
  )

set(_mumps_required_vars)
foreach(_precision ${MUMPS_PRECISIONS})
  string(TOUPPER "${_precision}" _u_precision)
  find_library(MUMPS_LIBRARY_${_u_precision}MUMPS NAMES ${_precision}mumps${MUMPS_PREFIX}
    PATHS "${MUMPS_DIR}"
    ENV MUMPS_DIR
    PATH_SUFFIXES lib
    )
  mark_as_advanced(MUMPS_LIBRARY_${_u_precision}MUMPS)
  list(APPEND _mumps_required_vars MUMPS_LIBRARY_${_u_precision}MUMPS)

  list(APPEND MUMPS_LIBRARIES_ALL ${MUMPS_LIBRARY_${_u_precision}MUMPS})
endforeach()

#===============================================================================
include(FindPackageHandleStandardArgs)
if(CMAKE_VERSION VERSION_GREATER 2.8.12)
  if(MUMPS_INCLUDE_DIR)
    file(STRINGS ${MUMPS_INCLUDE_DIR}/dmumps_c.h _versions
      REGEX "^#define MUMPS_VERSION .*")
    foreach(_ver ${_versions})
      string(REGEX MATCH "MUMPS_VERSION *\"([0-9.]+)\"" _tmp "${_ver}")
      set(_mumps_VERSION ${CMAKE_MATCH_1})
    endforeach()
    set(MUMPS_VERSION "${_mumps_VERSION}" CACHE INTERNAL "")
  endif()

  find_package_handle_standard_args(Mumps
    REQUIRED_VARS
      ${_mumps_required_vars}
      MUMPS_INCLUDE_DIR
    VERSION_VAR
      MUMPS_VERSION)
else()
  find_package_handle_standard_args(Mumps DEFAULT_MSG
    ${_mumps_required_vars} MUMPS_INCLUDE_DIR)
endif()

#===============================================================================
set(_mumps_dep_symbol_BLAS ${_first_precision}gemm)
set(_mumps_dep_symbol_ScaLAPACK numroc)
set(_mumps_dep_symbol_MPI mpi_send)
set(_mumps_dep_symbol_Scotch scotchfstratinit)
set(_mumps_dep_symbol_Scotch_ptscotch scotchfdgraphexit)
set(_mumps_dep_symbol_Scotch_esmumps esmumps)
set(_mumps_dep_symbol_mumps_common mumps_abort)
set(_mumps_dep_symbol_pord SPACE_ordering)
set(_mumps_dep_symbol_METIS metis_nodend)
set(_mumps_dep_symbol_ParMETIS ParMETIS_V3_NodeND)

 # added for fucking macosx that cannot fail at link
set(_mumps_run_dep_symbol_mumps_common mumps_fac_descband)
set(_mumps_run_dep_symbol_MPI mpi_abort)

set(_mumps_dep_link_BLAS BLAS_LIBRARIES)
set(_mumps_dep_link_ScaLAPACK SCALAPACK_LIBRARIES)
set(_mumps_dep_link_MPI MPI_Fortran_LIBRARIES)
set(_mumps_dep_link_Scotch SCOTCH_LIBRARIES)
set(_mumps_dep_link_Scotch_ptscotch SCOTCH_LIBRARY_PTSCOTCH)
set(_mumps_dep_link_Scotch_esmumps SCOTCH_LIBRARY_ESMUMPS)
set(_mumps_dep_link_mumps_common MUMPS_LIBRARY_COMMON)
set(_mumps_dep_link_pord MUMPS_LIBRARY_PORD)
set(_mumps_dep_link_METIS METIS_LIBRARY)
set(_mumps_dep_link_ParMETIS PARMETIS_LIBRARY)

set(_mumps_dep_comp_Scotch_ptscotch COMPONENTS ptscotch)
set(_mumps_dep_comp_Scotch_ptscotch COMPONENTS esmumps)

set(_mumps_potential_dependencies mumps_common pord BLAS ScaLAPACK MPI
  Scotch Scotch_ptscotch Scotch_esmumps METIS ParMETIS)
#===============================================================================

if(MUMPS_LIBRARY_${_u_first_precision}MUMPS MATCHES ".*${_first_precision}mumps.*${CMAKE_STATIC_LIBRARY_SUFFIX}")
  # Assuming mumps was compiled as a static library
  set(MUMPS_LIBRARY_TYPE STATIC CACHE INTERNAL "" FORCE)

  if (CMAKE_Fortran_COMPILER MATCHES ".*gfortran")
    set(_compiler_specific gfortran)
  elseif (CMAKE_Fortran_COMPILER MATCHES ".*ifort")
    set(_compiler_specific ifcore)
  else()
    message("Compiler ${CMAKE_Fortran_COMPILER} is not known, you will probably "
      "have to add semething instead of this message to be able to test mumps "
      "install")
  endif()
else()
  set(MUMPS_LIBRARY_TYPE SHARED CACHE INTERNAL "" FORCE)
endif()

try_run(_mumps_run _mumps_compiles "${_mumps_test_dir}" "${_mumps_test_dir}/mumps_test_code.c"
  CMAKE_FLAGS "-DINCLUDE_DIRECTORIES:STRING=${MUMPS_INCLUDE_DIR}"
  LINK_LIBRARIES ${MUMPS_LIBRARIES_ALL} ${_compiler_specific}
  RUN_OUTPUT_VARIABLE _run
  COMPILE_OUTPUT_VARIABLE _out)

foreach(_pdep ${_mumps_potential_dependencies})
  if(_mumps_compiles AND NOT _mumps_run STREQUAL "FAILED_TO_RUN")
    break()
  endif()
  
  if(_out MATCHES "${_mumps_dep_symbol_${_pdep}}" OR (DEFINED _mumps_run_dep_symbol_${_pdep} AND _run MATCHES "${_mumps_run_dep_symbol_${_pdep}}"))
    if(_pdep STREQUAL "mumps_common")
      find_library(MUMPS_LIBRARY_COMMON mumps_common${MUMPS_PREFIX}
        PATHS "${MUMPS_DIR}"
        ENV MUMPS_DIR
        PATH_SUFFIXES lib
        )

      if(NOT TARGET MUMPS::common)
        add_library(MUMPS::common ${MUMPS_LIBRARY_TYPE} IMPORTED GLOBAL)
      endif()
      set_target_properties(MUMPS::common PROPERTIES
        IMPORTED_LOCATION                 "${MUMPS_LIBRARY_COMMON}"
        INTERFACE_INCLUDE_DIRECTORIES     "${MUMPS_INCLUDE_DIR}"
        IMPORTED_LINK_INTERFACE_LANGUAGES "C;Fortran")
      list(APPEND _mumps_interface_link MUMPS::common)
    elseif(_pdep STREQUAL "pord")
      find_library(MUMPS_LIBRARY_PORD pord${MUMPS_PREFIX}
        PATHS "${MUMPS_DIR}"
        ENV MUMPS_DIR
        PATH_SUFFIXES lib
        )
      if(NOT TARGET MUMPS::pord)
        add_library(MUMPS::pord   ${MUMPS_LIBRARY_TYPE} IMPORTED GLOBAL)
      endif()
      #TODO adapt it for windows and dlls (check FindGSL as an example)
      set_target_properties(MUMPS::pord PROPERTIES
        IMPORTED_LOCATION                 "${MUMPS_LIBRARY_PORD}"
        INTERFACE_INCLUDE_DIRECTORIES     "${MUMPS_INCLUDE_DIR}"
        IMPORTED_LINK_INTERFACE_LANGUAGES "C")
      list(APPEND _mumps_interface_link MUMPS::pord)
    elseif(_pdep MATCHES "Scotch")
      find_package(Scotch REQUIRED ${_mumps_dep_comp_${_pdep}})
      list(APPEND _mumps_interface_link ${${_mumps_dep_link_${_pdep}}})
    else()
      find_package(${_pdep} REQUIRED)
      list(APPEND _mumps_interface_link ${${_mumps_dep_link_${_pdep}}})
    endif()

    list(APPEND MUMPS_LIBRARIES_ALL ${${_mumps_dep_link_${_pdep}}})

    try_run(_mumps_run _mumps_compiles "${_mumps_test_dir}" "${_mumps_test_dir}/mumps_test_code.c"
      CMAKE_FLAGS "-DINCLUDE_DIRECTORIES:STRING=${MUMPS_INCLUDE_DIR}"
      LINK_LIBRARIES ${MUMPS_LIBRARIES_ALL} ${_compiler_specific}
      RUN_OUTPUT_VARIABLE _run
      COMPILE_OUTPUT_VARIABLE _out)
  endif()
endforeach()

foreach(_precision ${MUMPS_PRECISIONS})
  string(TOUPPER "${_precision}" _u_precision)
  set(_target MUMPS::${_precision}mumps)
  if(NOT TARGET ${_target})
    add_library(${_target} ${MUMPS_LIBRARY_TYPE} IMPORTED GLOBAL)
  endif()
  set_target_properties(${_target} PROPERTIES
    IMPORTED_LOCATION                 "${MUMPS_LIBRARY_${_u_precision}MUMPS}"
    INTERFACE_INCLUDE_DIRECTORIES     "${MUMPS_INCLUDE_DIR}"
    IMPORTED_LINK_INTERFACE_LANGUAGES "C;Fortran"
    INTERFACE_LINK_LIBRARIES          "${_mumps_interface_link}")
endforeach()

set(MUMPS_LIBRARIES ${MUMPS_LIBRARIES_ALL} CACHE INTERNAL "" FORCE)
