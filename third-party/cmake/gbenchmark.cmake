set(_working_dir ${PROJECT_BINARY_DIR}/third-party/src/gbenchmark-download)

configure_file(${PROJECT_SOURCE_DIR}/third-party/gbenchmark.cmake.in ${_working_dir}/CMakeLists.txt)

set(GBENCHMARK_ROOT ${PROJECT_BINARY_DIR}/third-party)
find_package(benchmark QUIET)

if(NOT BENCHMARK_FOUND)
  if(NOT EXISTS ${PROJECT_SOURCE_DIR}/third-party/google-benchmark)
    message(STATUS "Downloading google-benchmark")
    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
      RESULT_VARIABLE result WORKING_DIRECTORY ${_working_dir}
      OUTPUT_FILE ${_working_dir}/configure-out.log
      ERROR_FILE ${_working_dir}/configure-error.log)

    if(result)
      message(SEND_ERROR "CMake step for google-benchmark failed: ${result}")
      return()
    endif()

    execute_process(COMMAND ${CMAKE_COMMAND} --build .
      RESULT_VARIABLE result WORKING_DIRECTORY ${_working_dir}
      OUTPUT_FILE ${_working_dir}/build-out.log
      ERROR_FILE ${_working_dir}/build-error.log)

    if(result)
      message(SEND_ERROR "Downloading google-benchmark failed: ${result}")
      return()
    endif()
  endif()

  set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE INTERNAL "" FORCE)
  set(BENCHMARK_ENABLE_TESTING OFF CACHE INTERNAL "" FORCE)
  add_subdirectory(third-party/google-benchmark)

  set_property(TARGET benchmark PROPERTY CXX_STANDARD 14)
  add_library(benchmark::benchmark ALIAS benchmark)

  set(gbenchmark_FOUND TRUE CACHE INTERNAL "" FORCE)
  set(GBENCHMARK_INCLUDE_DIRS third-party/google-benchmark/googletest/include CACHE INTERNAL "" FORCE)
  set(GBENCHMARK_LIBRARIES benchmark::benchmark CACHE INTERNAL "" FORCE)

  mask_package_options(BENCHMARK)
  mark_as_advanced(gbenchmark_DIR LIBRT LLVM_FILECHECK_EXE)
endif()
