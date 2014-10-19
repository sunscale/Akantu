#===============================================================================
# @file   90_scotch.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Nov 21 2011
# @date last modification: Thu Jul 10 2014
#
# @brief  package description for scotch
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

set(AKANTU_SCOTCH_FILES
  mesh_utils/mesh_partition/mesh_partition_scotch.cc
  )

if(AKANTU_SCOTCH_ON OR AKANTU_PTSCOTCH_ON)
  set(AKANTU_PARTITIONER_ON ON)
else()
  set(AKANTU_PARTITIONER_ON OFF)
endif()

set(AKANTU_SCOTCH_TESTS
  test_mesh_partitionate_scotch
  test_mesh_partitionate_scotch_advanced
  )


option(AKANTU_USE_SCOTCH "Add Scotch support in akantu" OFF)
option(AKANTU_USE_THIRD_PARTY_SCOTCH "Use the third-party Scotch instead of the one from the system" OFF)
mark_as_advanced(AKANTU_USE_THIRD_PARTY_SCOTCH)

set(SCOTCH_VERSION "5.1.12b")
set(SCOTCH_ARCHIVE_HASH "MD5=e13b49be804755470b159d7052764dc0")
set(SCOTCH_ARCHIVE ${PROJECT_SOURCE_DIR}/third-party/scotch_${SCOTCH_VERSION}_esmumps.tar.gz)
if(NOT EXISTS ${SCOTCH_ARCHIVE})
  set(SCOTCH_ARCHIVE https://gforge.inria.fr/frs/download.php/28978/scotch_${SCOTCH_VERSION}_esmumps.tar.gz)
endif()


if(AKANTU_USE_THIRD_PARTY_SCOTCH AND AKANTU_USE_SCOTCH)
  if(TARGET Scotch)
    return()
  endif()

  find_package(BISON)
  find_package(FLEX)
  find_package(ZLIB)

  if (AKANTU_USE_OBSOLETE_GETTIMEOFDAY)
    set (SCOTCH_TIMMING_OPTION -DCOMMON_TIMING_OLD)
  endif()

  if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(SCOTCH_ARCHITECTURE -DIDXSIZE64)
  else()
    set(SCOTCH_ARCHITECTURE)
  endif()

  set(SCOTCH_DIR ${PROJECT_BINARY_DIR}/third-party)
  configure_file(${PROJECT_SOURCE_DIR}/third-party/scotch_${SCOTCH_VERSION}_make.inc.cmake
    ${SCOTCH_DIR}/scotch_make.inc)

  include(ExternalProject)

  ExternalProject_Add(Scotch
    PREFIX ${SCOTCH_DIR}
    URL ${SCOTCH_ARCHIVE}
    URL_HASH ${SCOTCH_ARCHIVE_HASH}
    TLS_VERIFY FALSE
    PATCH_COMMAND patch -p1 < ${PROJECT_SOURCE_DIR}/third-party/scotch_${SCOTCH_VERSION}.patch
    CONFIGURE_COMMAND cmake -E copy ${SCOTCH_DIR}/scotch_make.inc src/Makefile.inc
    BUILD_IN_SOURCE 1
    BUILD_COMMAND make -C src
    INSTALL_COMMAND prefix=<INSTALL_DIR> make -C src install
    )

  set_third_party_shared_libirary_name(SCOTCH_LIBRARY scotch)
  set_third_party_shared_libirary_name(SCOTCH_LIBRARY_ERR     scotcherr)
  set_third_party_shared_libirary_name(SCOTCH_LIBRARY_ERREXIT scotcherrexit)
  set_third_party_shared_libirary_name(SCOTCH_LIBRARY_ESMUMPS esmumps)

  set(SCOTCH_INCLUDE_DIR     ${SCOTCH_DIR}/include CACHE PATH "" FORCE)

  mark_as_advanced(SCOTCH_LIBRARY SCOTCH_LIBRARY_ERR SCOTCH_LIBRARY_ERREXIT SCOTCH_LIBRARY_ESMUMPS
    SCOTCH_INCLUDE_DIR)

  set(SCOTCH_LIBRARIES_ALL ${SCOTCH_LIBRARY} ${SCOTCH_LIBRARY_ERR})
  set(SCOTCH_LIBRARIES ${SCOTCH_LIBRARIES_ALL} CACHE INTERNAL "Libraries for scotch" FORCE)

  list(APPEND AKANTU_EXTERNAL_LIBRARIES ${SCOTCH_LIBRARIES})
  list(APPEND AKANTU_EXTERNAL_LIB_INCLUDE_DIR ${SCOTCH_INCLUDE_DIR})

  set(AKANTU_SCOTCH_INCLUDE_DIR ${SCOTCH_INCLUDE_DIR})
  set(AKANTU_SCOTCH_LIBRARIES ${SCOTCH_LIBRARIES})

  list(APPEND AKANTU_OPTION_LIST SCOTCH)
  set(SCOTCH_FOUND TRUE CACHE INTERNAL "" FORCE)
  set(AKANTU_SCOTCH ON)

  list(APPEND AKANTU_EXTRA_TARGET_DEPENDENCIES Scotch)
else()
  add_optional_external_package(Scotch "Add Scotch support in akantu" OFF)

  if(SCOTCH_INCLUDE_DIR)
    file(STRINGS ${SCOTCH_INCLUDE_DIR}/scotch.h SCOTCH_INCLUDE_CONTENT)
    string(REGEX MATCH "_cplusplus" _match ${SCOTCH_INCLUDE_CONTENT})
    if(_match)
      set(AKANTU_SCOTCH_NO_EXTERN ON)
      list(APPEND AKANTU_DEFINITIONS AKANTU_SCOTCH_NO_EXTERN)
    else()
      set(AKANTU_SCOTCH_NO_EXTERN OFF)
    endif()
  endif()
endif()

set(AKANTU_SCOTCH_DOCUMENTATION "
This package enables the use the \\href{http://www.labri.fr/perso/pelegrin/scotch/}{Scotch} library in
order to perform a graph partitioning leading to the domain
decomposition used within \\akantu

Under Ubuntu (14.04 LTS) the installation can be performed using the commands:
\\begin{command}
  > sudo apt-get install libscotch-dev
\\end{command}

If you activate the advanced option AKANTU\\_USE\\_THIRD\\_PARTY\\_SCOTCH the make system of akantu can automatically compile Scotch.

If the automated download fails due to a SSL access not supported by your version of CMake please download the file \\href{https://gforge.inria.fr/frs/download.php/28978/scotch\\_${SCOTCH_VERSION}\\_esmumps.tar.gz}{scotch\\_${SCOTCH_VERSION}\\_esmumps.tar.gz}  and then place it in the directory \\shellcode{<akantu source>/third-party}
" )
