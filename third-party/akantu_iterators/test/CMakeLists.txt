enable_testing()

add_executable(test_iterators test_iterators.cc)
target_link_libraries(test_iterators
  PRIVATE akantu_iterators GTest::GTest)

include(GoogleTest)
gtest_discover_tests(test_iterators)
