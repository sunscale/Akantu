#===============================================================================
# @file   scotch.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Mar 30 2015
# @date last modification: Wed Nov 11 2015
#
# @brief  build script for scotch
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

if(TARGET Scotch)
  return()
endif()

if(NOT EXISTS ${PROJECT_SOURCE_DIR}/third-party/${SCOTCH_ARCHIVE})
  set(_scotch_download_command
    URL ${SCOTCH_URL}
#    URL_HASH ${SCOTCH_ARCHIVE_HASH}
    TLS_VERIFY FALSE
    )
else()
  set(_scotch_download_command
    URL ${PROJECT_SOURCE_DIR}/third-party/${SCOTCH_ARCHIVE}
    URL_HASH ${SCOTCH_ARCHIVE_HASH})
endif()

if(CMAKE_VERSION VERSION_GREATER 3.1)
  set(_extra_options
    DOWNLOAD_NO_PROGRESS 1
    EXCLUDE_FROM_ALL 1
    )
endif()

find_package(BISON REQUIRED)
find_package(FLEX REQUIRED)
find_package(ZLIB)

#if(ZLIB_FOUND)
#  set(_zlib_cflags "-DCOMMON_FILE_COMPRESS_GZ -I${ZLIB_INCLUDE_DIR}")
#  set(_zlib_ldflags "${ZLIB_LIBRARY}")
#endif()

if(NOT CMAKE_SYSTEM_NAME STREQUAL "Windows")
  if (AKANTU_USE_OBSOLETE_GETTIMEOFDAY)
    set(_timing_cflags -DCOMMON_TIMING_OLD)
  endif()
  set(_system_cflags "-DCOMMON_PTHREAD -DSCOTCH_PTHREAD ${_timing_cflags}")
  set(_system_ldflags "-lpthread")
else()
  set(_system_cflags "-DCOMMON_RANDOM_RAND -DCOMMON_WINDOWS -DCOMMON_STUB_FORK -D'pipe(pfds)=_pipe(pfds,1024,0x8000)'")
endif()

if(CMAKE_SIZEOF_VOID_P EQUAL 8)
  set(_architecture_cflags -DIDXSIZE64)
endif()

math(EXPR _n "${AKANTU_INTEGER_SIZE} * 8")
if(NOT _n EQUAL 32)
  set(_num_size_cflags "-DINTSIZE${_n}")
endif()

if(HAVE_STDINT_H)
  set(_stdint -DHAVE_STDINT_H)
endif()

set(AKANTU_SCOTCH_CFLAGS "-O3 -w -fPIC -Drestrict=__restrict -DCOMMON_RANDOM_FIXED_SEED -DSCOTCH_RENAME -DSCOTCH_RENAME_PARSER ${_zlib_cflags} ${_system_cflags} ${_architecture_cflags} ${_num_size_cflags} ${_stdint}")
set(AKANTU_SCOTCH_LDFLAGS "${_zlib_ldflags} ${_system_ldflags} -lm")

set(SCOTCH_DIR ${PROJECT_BINARY_DIR}/third-party)
configure_file(
  ${PROJECT_SOURCE_DIR}/third-party/scotch_${SCOTCH_VERSION}_make.inc.cmake
  ${SCOTCH_DIR}/scotch_make.inc)

include(ExternalProject)

ExternalProject_Add(Scotch
  PREFIX ${SCOTCH_DIR}
  ${_scotch_download_command}
  ${_extra_options}
  PATCH_COMMAND ${PATCH_COMMAND} -p1 < ${PROJECT_SOURCE_DIR}/third-party/scotch_${SCOTCH_VERSION}.patch
  CONFIGURE_COMMAND ${CMAKE_COMMAND} -E copy ${SCOTCH_DIR}/scotch_make.inc src/Makefile.inc
  BUILD_IN_SOURCE 1
  BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} -C src
  INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} prefix=<INSTALL_DIR> -C src install
  LOG_DOWNLOAD 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
  )

set_third_party_shared_libirary_name(SCOTCH_LIBRARY scotch)
set_third_party_shared_libirary_name(SCOTCH_LIBRARY_ERR     scotcherr)
set_third_party_shared_libirary_name(SCOTCH_LIBRARY_ERREXIT scotcherrexit)
set_third_party_shared_libirary_name(SCOTCH_LIBRARY_ESMUMPS esmumps)

set(SCOTCH_INCLUDE_DIR ${SCOTCH_DIR}/include CACHE PATH "" FORCE)

mark_as_advanced(
  SCOTCH_LIBRARY
  SCOTCH_LIBRARY_ERR
  SCOTCH_LIBRARY_ERREXIT
  SCOTCH_LIBRARY_ESMUMPS
  SCOTCH_INCLUDE_DIR)

set(SCOTCH_LIBRARIES_ALL ${SCOTCH_LIBRARY} ${SCOTCH_LIBRARY_ERR})
set(SCOTCH_LIBRARIES ${SCOTCH_LIBRARIES_ALL} CACHE INTERNAL "Libraries for scotch" FORCE)

package_add_extra_dependency(Scotch Scotch)
