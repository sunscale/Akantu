set(BLAS_DIR ${PROJECT_BINARY_DIR}/third-party)

configure_file(
  ${PROJECT_SOURCE_DIR}/third-party/blas_${BLAS_VERSION}_make.inc.cmake
  ${BLAS_DIR}/blas_make.inc @ONLY)

file(MAKE_DIRECTORY ${BLAS_DIR}/lib)

ExternalProject_Add(netlib-blas
  PREFIX ${BLAS_DIR}
  URL ${BLAS_ARCHIVE}
  CONFIGURE_COMMAND cmake -E copy ${BLAS_DIR}/blas_make.inc make.inc
  BUILD_IN_SOURCE 1
  BUILD_COMMAND ${CMAKE_MAKE_PROGRAM}
  INSTALL_COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SHARED_LIBRARY_PREFIX}blas${CMAKE_SHARED_LIBRARY_SUFFIX} <INSTALL_DIR>/lib
  LOG_DOWNLOAD 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
  )

package_add_extra_dependency(BLAS netlib-blas)

set_third_party_shared_libirary_name(BLAS_LIBRARIES blas)
