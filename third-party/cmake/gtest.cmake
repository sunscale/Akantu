set(_working_dir ${PROJECT_BINARY_DIR}/third-party/src/gtest-download)

configure_file(${PROJECT_SOURCE_DIR}/third-party/gtest.cmake.in ${_working_dir}/CMakeLists.txt)

set(GTEST_ROOT ${PROJECT_BINARY_DIR}/third-party)
find_package(GTest QUIET)


if(NOT GTEST_FOUND)
  message(STATUS "Downloading googletest")
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result WORKING_DIRECTORY ${_working_dir}
    OUTPUT_FILE ${_working_dir}/configure-out.log
    ERROR_FILE ${_working_dir}/configure-error.log)

  if(result)
    message(SEND_ERROR "CMake step for googletest failed: ${result}")
    return()
  endif()

  message(STATUS "Building googletest")
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result WORKING_DIRECTORY ${_working_dir}
    OUTPUT_FILE ${_working_dir}/build-out.log
    ERROR_FILE ${_working_dir}/build-error.log)

  if(result)
    message(SEND_ERROR "Build step for googletest failed: ${result}")
    return()
  endif()
endif()

set(GTEST_ROOT ${PROJECT_BINARY_DIR}/third-party)
find_package(GTest REQUIRED)
