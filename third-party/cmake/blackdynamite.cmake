if(${PROJECT_SOURCE_DIR}/third-party/${BLACKDYNAMITE_ARCHIVE})
  set(_blackdynamite_download_command
    URL ${PROJECT_SOURCE_DIR}/third-party/${BLACKDYNAMITE_ARCHIVE})
else()
  set(_blackdynamite_download_command
    GIT_REPOSITORY ${BLACKDYNAMITE_GIT}
    GIT_TAG ${BLACKDYNAMITE_VERSION}
    )
endif()

if(CMAKE_VERSION VERSION_GREATER 3.1)
  set(_extra_options 
    UPDATE_DISCONNECTED 1
    DOWNLOAD_NO_PROGRESS 1
    EXCLUDE_FROM_ALL 1
    )
endif()

ExternalProject_Add(blackdynamite
  PREFIX ${PROJECT_BINARY_DIR}/third-party
  ${_blackdynamite_download_command}
  ${_extra_options}
  CMAKE_ARGS <SOURCE_DIR>/
  CMAKE_CACHE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE} -DCMAKE_CXX_COMPILER:PATH=${CMAKE_CXX_COMPILER}
  BUILD_COMMAND make
  INSTALL_COMMAND make install
  LOG_DOWNLOAD 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
  )


set_third_party_shared_libirary_name(BLACKDYNAMITE_LIBRARIES blackdynamite)
set(BLACKDYNAMITE_INCLUDE_DIR "${PROJECT_BINARY_DIR}/third-party/include/blackdynamite" CACHE PATH "")
mark_as_advanced(
  BLACKDYNAMITE_LIBRARIES
  BLACKDYNAMITE_INCLUDE_DIR
  )

package_add_extra_dependency(BlackDynamite blackdynamite)
