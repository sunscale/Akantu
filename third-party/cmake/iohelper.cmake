if(NOT EXISTS ${PROJECT_SOURCE_DIR}/third-party/${IOHELPER_ARCHIVE})
  set(_iohelper_download_command
    GIT_REPOSITORY ${IOHELPER_GIT}
    GIT_TAG ${IOHELPER_VERSION}
    )
else()
  set(_iohelper_download_command
    URL ${PROJECT_SOURCE_DIR}/third-party/${IOHELPER_ARCHIVE}
    )
endif()

if(CMAKE_VERSION VERSION_GREATER 3.1)
  set(_extra_options 
    UPDATE_DISCONNECTED 1
    DOWNLOAD_NO_PROGRESS 1
    EXCLUDE_FROM_ALL 1
    )
endif()

ExternalProject_Add(IOHelper
  PREFIX ${PROJECT_BINARY_DIR}/third-party
  ${_iohelper_download_command}
  ${_extra_options}
  CMAKE_ARGS <SOURCE_DIR>/
  CMAKE_CACHE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE} -DCMAKE_CXX_COMPILER:PATH=${CMAKE_CXX_COMPILER}
  LOG_DOWNLOAD 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
  )

set_third_party_shared_libirary_name(IOHELPER_LIBRARIES iohelper)
set(IOHELPER_INCLUDE_DIR "${PROJECT_BINARY_DIR}/third-party/include/iohelper" CACHE PATH "IOHelper include directory")

mark_as_advanced(
  IOHELPER_LIBRARIES
  IOHELPER_INCLUDE_DIR
  )

package_add_extra_dependency(IOHelper IOHelper)
