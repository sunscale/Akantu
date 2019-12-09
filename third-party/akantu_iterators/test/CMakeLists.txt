enable_testing()

add_executable(test_iterators test_iterators.cc)
target_link_libraries(test_iterators
  PRIVATE akantu_iterators GTest::GTest)

if ((CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "5.2") OR
    CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  target_compile_options(test_iterators PRIVATE -fsanitize=address -fsanitize=leak -fsanitize=undefined -fno-omit-frame-pointer)
  target_link_options(test_iterators PRIVATE -fsanitize=address -fsanitize=leak -fsanitize=undefined -fno-omit-frame-pointer)
endif()

include(GoogleTest)
gtest_discover_tests(test_iterators EXTRA_ARGS "--gtest_output=xml")
