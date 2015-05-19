ExternalProject_Add(IOHelper
  PREFIX ${PROJECT_BINARY_DIR}/third-party
  GIT_REPOSITORY ${IOHELPER_GIT}
  CMAKE_ARGS <SOURCE_DIR>/
  CMAKE_CACHE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE} -DCMAKE_CXX_COMPILER:PATH=${CMAKE_CXX_COMPILER}
  TIMEOUT 6
  )

set_third_party_shared_libirary_name(IOHELPER_LIBRARIES iohelper)
set(IOHELPER_INCLUDE_DIR "${PROJECT_BINARY_DIR}/third-party/include" CACHE PATH "IOHelper include directory")

mark_as_advanced(
  IOHELPER_LIBRARIES
  IOHELPER_INCLUDE_DIR
  )

package_set_libraries(IOHelper ${IOHELPER_LIBRARIES})
package_set_include_dir(IOHelper ${PROJECT_BINARY_DIR}/third-party/include/iohelper)

package_add_extra_dependency(IOHelper IOHelper)
