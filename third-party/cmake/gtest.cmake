set(_working_dir ${PROJECT_BINARY_DIR}/third-party/src/gtest-download)

configure_file(${PROJECT_SOURCE_DIR}/third-party/gtest.cmake.in ${_working_dir}/CMakeLists.txt)

set(GTEST_ROOT ${PROJECT_BINARY_DIR}/third-party)

if(NOT EXISTS ${PROJECT_SOURCE_DIR}/third-party/google-test)
  message(STATUS "Downloading googletest")
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result WORKING_DIRECTORY ${_working_dir}
    OUTPUT_FILE ${_working_dir}/configure-out.log
    ERROR_FILE ${_working_dir}/configure-error.log)

  if(result)
    message(SEND_ERROR "CMake step for googletest failed: ${result}")
    return()
  endif()

  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result WORKING_DIRECTORY ${_working_dir}
    OUTPUT_FILE ${_working_dir}/build-out.log
    ERROR_FILE ${_working_dir}/build-error.log)

  if(result)
    message(SEND_ERROR "Downloading googletest failed: ${result}")
    return()
  endif()
endif()

set(gtest_build_tests OFF CACHE INTERNAL "" FORCE)
set(BUILD_GTEST ON CACHE INTERNAL "" FORCE)
set(BUILD_GMOCK OFF CACHE INTERNAL "" FORCE)

set(Python_ADDITIONAL_VERSIONS ${AKANTU_PREFERRED_PYTHON_VERSION})
add_subdirectory(third-party/google-test)
set_property(TARGET gtest_main PROPERTY CXX_STANDARD 14)
set_property(TARGET gtest PROPERTY CXX_STANDARD 14)

add_library(GTest::Main ALIAS gtest_main)
add_library(GTest::GTest ALIAS gtest)

set(gtest_FOUND TRUE CACHE INTERNAL "" FORCE)
set(GTEST_INCLUDE_DIRS third-party/google-test/googletest/include CACHE INTERNAL "" FORCE)
set(GTEST_LIBRARIES GTest::Main CACHE INTERNAL "" FORCE)

mark_as_advanced(INSTALL_GTEST)
mask_package_options(gtest)
