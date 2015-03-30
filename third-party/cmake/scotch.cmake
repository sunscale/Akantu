set(SCOTCH_ARCHIVE ${PROJECT_SOURCE_DIR}/third-party/scotch_${SCOTCH_VERSION}_esmumps.tar.gz)
if(NOT EXISTS ${SCOTCH_ARCHIVE})
  set(SCOTCH_ARCHIVE ${SCOTCH_URL})
endif()

if(TARGET Scotch)
  return()
endif()

find_package(BISON REQUIRED)
find_package(FLEX REQUIRED)
find_package(ZLIB REQUIRED)

if (AKANTU_USE_OBSOLETE_GETTIMEOFDAY)
  set (SCOTCH_TIMMING_OPTION -DCOMMON_TIMING_OLD)
endif()

if(CMAKE_SIZEOF_VOID_P EQUAL 8)
  set(SCOTCH_ARCHITECTURE -DIDXSIZE64)
else()
  set(SCOTCH_ARCHITECTURE)
endif()

set(SCOTCH_DIR ${PROJECT_BINARY_DIR}/third-party)
configure_file(
  ${PROJECT_SOURCE_DIR}/third-party/scotch_${SCOTCH_VERSION}_make.inc.cmake
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
  BUILD_COMMAND MPICH_CC=${CMAKE_C_COMPILER} ${CMAKE_MAKE_PROGRAM} -C src
  INSTALL_COMMAND prefix=<INSTALL_DIR> ${CMAKE_MAKE_PROGRAM} -C src install
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

package_set_libraries(Scotch   ${SCOTCH_LIBRARIES})
package_set_include_dir(Scotch ${SCOTCH_INCLUDE_DIR})

package_add_extra_dependency(Scotch Scotch)