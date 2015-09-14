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
find_package(ZLIB REQUIRED)

if (AKANTU_USE_OBSOLETE_GETTIMEOFDAY)
  set (SCOTCH_TIMMING_OPTION -DCOMMON_TIMING_OLD)
endif()

if(CMAKE_SIZEOF_VOID_P EQUAL 8)
  set(SCOTCH_ARCHITECTURE -DIDXSIZE64)
else()
  set(SCOTCH_ARCHITECTURE)
endif()

math(EXPR _n "${AKANTU_INTEGER_SIZE} * 8")
set(SCOTCH_NUM_SIZE "-DINTSIZE${_n}")

set(SCOTCH_DIR ${PROJECT_BINARY_DIR}/third-party)
configure_file(
  ${PROJECT_SOURCE_DIR}/third-party/scotch_${SCOTCH_VERSION}_make.inc.cmake
  ${SCOTCH_DIR}/scotch_make.inc)

include(ExternalProject)

ExternalProject_Add(Scotch
  PREFIX ${SCOTCH_DIR}
  ${_scotch_download_command}
  ${_extra_options}
  PATCH_COMMAND patch -p1 < ${PROJECT_SOURCE_DIR}/third-party/scotch_${SCOTCH_VERSION}.patch
  CONFIGURE_COMMAND cmake -E copy ${SCOTCH_DIR}/scotch_make.inc src/Makefile.inc
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
