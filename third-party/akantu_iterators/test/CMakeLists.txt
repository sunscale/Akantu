enable_testing()

add_executable(test_akantu_iterators test_gtest_main.cc test_tuples.cc test_str.cc test_iterators.cc)
target_link_libraries(test_akantu_iterators
  PRIVATE akantu_iterators GTest::GTest)
target_compile_definitions(test_akantu_iterators PRIVATE -DAKANTU_ITERATORS_NAMESPACE=aka)

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  target_compile_options(test_akantu_iterators PRIVATE -Wall -Wextra -pedantic -Weffc++)
endif()

if ((CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "5.2") OR
    CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  target_compile_options(test_akantu_iterators PRIVATE -fsanitize=address -fsanitize=leak -fsanitize=undefined -fno-omit-frame-pointer)
  target_link_options(test_akantu_iterators PRIVATE -fsanitize=address -fsanitize=leak -fsanitize=undefined -fno-omit-frame-pointer)
endif()

include(GoogleTest)
gtest_discover_tests(test_akantu_iterators)
