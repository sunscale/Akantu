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
set(SCALAPACK_VENDOR "Auto" CACHE
  STRING "Vendor for scalapack (Auto, Netlib, Intel_(i?)lp64_(openmpi|intelmpi|sgimpt)")
mark_as_advanced(SCALAPACK_VENDOR)
set_property(CACHE SCALAPACK_VENDOR PROPERTY STRINGS
  Auto Netlib
  Intel_lp64_openmpi Intel_lp64_intelmpi Intel_lp64_sgimpt
  Intel_ilp64_openmpi Intel_ilp64_intelmpi Intel_ilp64_sgimpt)


macro(scalapack_find_library prefix target name list_libraries list_headers)
  foreach(_lib ${list_libraries})
    find_library(${prefix}_${_lib}_LIBRARY NAMES ${_lib}
      PATHS ${prefix}_DIR
      )
    mark_as_advanced(${prefix}_${_lib}_LIBRARY)

    if(${prefix}_${_lib}_LIBRARY)
      list(APPEND ${prefix}_libraries ${${prefix}_${_lib}_LIBRARY})

      get_filename_component(_ext ${${prefix}_${_lib}_LIBRARY} EXT)
      if(NOT TARGET ${target})
	if("${_ext}" STREQUAL "${CMAKE_SHARED_LIBRARY_SUFFIX}")
	  add_library(${target} SHARED IMPORTED)
	  get_filename_component(_soname ${${prefix}_${_lib}_LIBRARY} NAME)
	  set_property(TARGET ${target} PROPERTY IMPORTED_SONAME ${_soname})
	else()
	  add_library(${target}_${name} STATIC IMPORTED)
	endif()
	set_property(TARGET ${target} PROPERTY
	    IMPORTED_LOCATION ${${prefix}_${_lib}_LIBRARY})
      else()
	if("${_ext}" STREQUAL "${CMAKE_SHARED_LIBRARY_SUFFIX}")
	  set_property(TARGET ${target} APPEND PROPERTY
	    IMPORTED_LINK_DEPENDENT_LIBRARIES ${${prefix}_${_lib}_LIBRARY}
	    )
	else()
	  set_property(TARGET ${target} APPEND PROPERTY
	    IMPORTED_LINK_INTERFACE_LIBRARIES ${${prefix}_${_lib}_LIBRARY}
	    )
	endif()
      endif()
    else()
      unset(${prefix}_${_lib}_LIBRARY CACHE)
    endif()
  endforeach()

  if(${prefix}_libraries)
    foreach(_hdr ${list_headers})
      get_filename_component(_hdr_name ${_hdr} NAME_WE)
      find_path(${prefix}_${_hdr_name}_INCLUDE_DIR NAMES ${_hdr}
	PATHS ${prefix}_DIR)
      mark_as_advanced(${prefix}_${_hdr_name}_INCLUDE_DIR)

      if(${prefix}_${_hdr_name}_INCLUDE_DIR)
	list(APPEND ${prefix}_include_dir ${${prefix}_${_hdr_name}_INCLUDE_DIR})
	set_property(TARGET ${target} APPEND PROPERTY
	  INTERFACE_INCLUDE_DIRECTORIES ${${prefix}_${_lib}_INCLUDE_DIR}
	  )
      else()
	unset(${prefix}_${_lib}_INCLUDE_DIR CACHE)
      endif()
    endforeach()
  endif()
endmacro()

set(SCALAPACK_libraries)
set(SCALAPACK_INCLUDE_DIR)

if(SCALAPACK_VENDOR STREQUAL "Auto" OR SCALAPACK_VENDOR STREQUAL "Netlib")
  if(NOT SCALAPACK_libraries)
    scalapack_find_library(
      SCALAPACK
      ScaLAPACK
      "netlib"
      "scalapack;blacsC;blacsF77;blacs"
      ""
      )
  endif()
endif()

foreach(_precision lp64 ilp64)
  foreach(_mpi intelmpi openmpi sgimpt)
    if(NOT SCALAPACK_libraries)
      if(SCALAPACK_VENDOR STREQUAL "Auto" OR SCALAPACK_VENDOR STREQUAL "Intel_${_precision}_${_mpi}")
	if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
	  set(_mkl_common "mkl_intel_${_precision}")

	else()
	  set(_mkl_common "mkl_gf_${_precision}")
	endif()
	scalapack_find_library(
	  SCALAPACK
	  ScaLAPACK
	  "intel_${_precision}_${_mpi}"
	  "mkl_scalapack_${_precision};${_mkl_common};mkl_sequential;mkl_core;mkl_blacs_${_mpi}_${_precision}"
	  "mkl_scalapack.h"
	  )

	if(SCALAPACK_libraries AND _precision STREQUAL "ilp64")
	  set_property(TARGET ${target} APPEND PROPERTY
	    INTERFACE_COMPILE_DEFINITIONS MKL_ILP64}
	    )
	endif()

	if(EXISTS ${SCALAPACK_include_dir}/mkl_version.h)
	  file(STRINGS ${SCALAPACK_include_dir}/mkl_version.h _versions
	    REGEX "^#define\ +__INTEL_MKL(_MINOR|_UPDATE)?__ .*")
	  foreach(_ver ${_versions})
	    string(REGEX MATCH "__INTEL_MKL(_MINOR|_UPDATE)?__ *([0-9.]+)" _tmp "${_ver}")
	    set(_mkl${CMAKE_MATCH_1} ${CMAKE_MATCH_2})
	  endforeach()
	  set(SCALAPACK_VERSION "mkl:${_mkl}.${_mkl_MINOR}.${_mkl_UPDATE}" CACHE INTERNAL "")
	endif()

      endif()
    endif()
  endforeach()
endforeach()

set(SCALAPACK_LIBRARIES ${SCALAPACK_libraries} CACHE INTERNAL "")
set(SCALAPACK_INCLUDE_DIR ${SCALAPACK_include_dir} CACHE INTERNAL "")

#===============================================================================
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ScaLAPACK
  REQUIRED_VARS SCALAPACK_LIBRARIES
  VERSION_VAR SCALAPACK_VERSION)
