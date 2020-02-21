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
set(MUMPS_PLAT)
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

find_path(MUMPS_INCLUDE_DIR ${_first_precision}mumps_c.h
  PATHS "${MUMPS_DIR}"
  ENV MUMPS_DIR
  PATH_SUFFIXES include
  )
mark_as_advanced(MUMPS_INCLUDE_DIR)

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


function(mumps_add_dependency _pdep _libs)
  string(TOUPPER ${_pdep} _u_pdep)
  if(_pdep STREQUAL "mumps_common")
    find_library(MUMPS_LIBRARY_COMMON mumps_common${MUMPS_PREFIX}
      PATHS "${MUMPS_DIR}"
      ENV MUMPS_DIR
      PATH_SUFFIXES lib
      )
    set(${_libs} ${MUMPS_LIBRARY_COMMON} PARENT_SCOPE)
    mark_as_advanced(MUMPS_LIBRARY_COMMON)
  elseif(_pdep STREQUAL "pord")
    find_library(MUMPS_LIBRARY_PORD pord${MUMPS_PREFIX}
      PATHS "${MUMPS_DIR}"
      ENV MUMPS_DIR
      PATH_SUFFIXES lib
      )
    set(${_libs} ${MUMPS_LIBRARY_PORD} PARENT_SCOPE)
    mark_as_advanced(MUMPS_LIBRARY_PORD)
  elseif(_pdep MATCHES "Scotch")
    find_package(Scotch REQUIRED ${ARGN} QUIET)
    if(ARGN)
      list(GET ARGN 1 _comp)
      string(TOUPPER ${_comp} _u_comp)
      set(${_libs} ${SCOTCH_LIBRARY_${_u_comp}} PARENT_SCOPE)
    else()
      set(${_libs} ${${_u_pdep}_LIBRARIES} PARENT_SCOPE)
    endif()
  elseif(_pdep MATCHES "MPI")
    if(MUMPS_PLAT STREQUAL "_seq")
      find_library(MUMPS_LIBRARY_MPISEQ mpiseq${MUMPS_PREFIX}
        PATHS "${MUMPS_DIR}"
        ENV MUMPS_DIR
        PATH_SUFFIXES lib
        )
      set(${_libs} ${MUMPS_LIBRARY_MPISEQ} PARENT_SCOPE)
      mark_as_advanced(MUMPS_LIBRARY_MPISEQ)
    else()
      find_package(MPI REQUIRED C Fortran QUIET)
      set(${_libs} ${MPI_C_LIBRARIES} ${MPI_Fortran_LIBRARIES} PARENT_SCOPE)
    endif()
  else()
    find_package(${_pdep} REQUIRED QUIET)
    set(${_libs} ${${_u_pdep}_LIBRARIES} ${${_u_pdep}_LIBRARY} PARENT_SCOPE)
  endif()
endfunction()

function(mumps_find_dependencies)
  set(_libraries_all ${MUMPS_LIBRARIES_ALL})
  set(_include_dirs ${MUMPS_INCLUDE_DIR})

  set(_mumps_test_dir "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}")
  file(READ ${CMAKE_CURRENT_LIST_DIR}/CheckFindMumps.c _output)
  file(WRITE "${_mumps_test_dir}/mumps_test_code.c"
    "#include <${_first_precision}mumps_c.h>
${_u_first_precision}MUMPS_STRUC_C id;

#define mumps_c ${_first_precision}mumps_c
#define Real ${_u_first_precision}MUMPS_REAL
")

  if(MUMPS_PLAT STREQUAL _seq)
    file(APPEND "${_mumps_test_dir}/mumps_test_code.c"
      "#define MUMPS_SEQ
")
  else()
    file(APPEND "${_mumps_test_dir}/mumps_test_code.c"
      "// #undef MUMPS_SEQ
")
        find_package(MPI REQUIRED)
    list(APPEND _compiler_specific ${MPI_C_LIBRARIES})
    list(APPEND _include_dirs ${MPI_C_INCLUDE_PATH} ${MPI_INCLUDE_DIR})
  endif()

  file(APPEND "${_mumps_test_dir}/mumps_test_code.c" "${_output}")

  #===============================================================================
  set(_mumps_dep_symbol_BLAS ${_first_precision}gemm)
  set(_mumps_dep_symbol_ScaLAPACK numroc)
  set(_mumps_dep_symbol_MPI mpi_send)
  set(_mumps_dep_symbol_Scotch SCOTCH_graphInit)
  set(_mumps_dep_symbol_Scotch_ptscotch scotchfdgraphexit)
  set(_mumps_dep_symbol_Scotch_esmumps esmumps)
  set(_mumps_dep_symbol_mumps_common mumps_abort)
  set(_mumps_dep_symbol_pord SPACE_ordering)
  set(_mumps_dep_symbol_METIS metis_nodend)
  set(_mumps_dep_symbol_ParMETIS ParMETIS_V3_NodeND)

  # added for fucking macosx that cannot fail at link
  set(_mumps_run_dep_symbol_mumps_common mumps_fac_descband)
  set(_mumps_run_dep_symbol_MPI mpi_bcast)
  set(_mumps_run_dep_symbol_ScaLAPACK idamax)

  set(_mumps_dep_comp_Scotch_ptscotch COMPONENTS ptscotch)
  set(_mumps_dep_comp_Scotch_esmumps COMPONENTS esmumps)

  set(_mumps_potential_dependencies mumps_common pord BLAS ScaLAPACK MPI
    Scotch Scotch_ptscotch Scotch_esmumps METIS ParMETIS)
  #===============================================================================

  set(_retry_try_run TRUE)
  set(_retry_count 0)

  # trying only as long as we add dependencies to avoid inifinte loop in case of an unkown dependency
  while (_retry_try_run AND _retry_count LESS 100)
    try_run(_mumps_run _mumps_compiles "${_mumps_test_dir}" "${_mumps_test_dir}/mumps_test_code.c"
      CMAKE_FLAGS "-DINCLUDE_DIRECTORIES:STRING=${_include_dirs}"
      LINK_LIBRARIES ${_libraries_all} ${_libraries_all} ${_compiler_specific}
      RUN_OUTPUT_VARIABLE _run
      COMPILE_OUTPUT_VARIABLE _out)

    set(_retry_compile FALSE)
    #message("COMPILATION outputs: \n${_out} \n RUN OUTPUT \n${_run}")
    if(_mumps_compiles AND NOT (_mumps_run STREQUAL "FAILED_TO_RUN"))
      break()
    endif()

    foreach(_pdep ${_mumps_potential_dependencies})
      #message("CHECKING ${_pdep}")
      set(_add_pdep FALSE)
      if (NOT _mumps_compiles AND
          _out MATCHES "undefined reference.*${_mumps_dep_symbol_${_pdep}}")
        set(_add_pdep TRUE)
        #message("NEED COMPILE ${_pdep}")
      elseif(_mumps_run STREQUAL "FAILED_TO_RUN" AND
          DEFINED _mumps_run_dep_symbol_${_pdep} AND
          _run MATCHES "${_mumps_run_dep_symbol_${_pdep}}")
        set(_add_pdep TRUE)
	#message("NEED RUN ${_pdep}")
      endif()

      if(_add_pdep)
        mumps_add_dependency(${_pdep} _libs ${_mumps_dep_comp_${_pdep}})
	#message("ADDING ${_libs}")
	if(NOT _libs)
	  message(FATAL_ERROR "MUMPS depends on ${_pdep} but no libraries where found")
	endif()
	list(APPEND _libraries_all ${_libs})
        set(_retry_try_run TRUE)
      endif()
    endforeach()

    math(EXPR _retry_count "${_retry_count} + 1")
  endwhile()

  if(_retry_count GREATER 10)
    message(FATAL_ERROR "Do not know what to do to link with mumps on your system, I give up!")
  endif()

  if(APPLE)
    # in doubt add some stuff because mumps was perhaps badly compiled
    mumps_add_dependency(pord _libs)
    list(APPEND _libraries_all ${_libs})
  endif()

  set(MUMPS_LIBRARIES_ALL ${_libraries_all} PARENT_SCOPE)
endfunction()

mumps_find_dependencies()

set(MUMPS_LIBRARIES ${MUMPS_LIBRARIES_ALL} CACHE INTERNAL "" FORCE)

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
    REQUIRED_VARS ${_mumps_required_vars}
                  MUMPS_INCLUDE_DIR
    VERSION_VAR MUMPS_VERSION
    )
else()
  find_package_handle_standard_args(Mumps DEFAULT_MSG
    ${_mumps_required_vars} MUMPS_INCLUDE_DIR)
endif()
