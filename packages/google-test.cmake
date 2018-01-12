package_declare(gtest EXTERNAL
  DESCRIPTION "Add GTest support for tests"
  SYSTEM AUTO third-party/cmake/gtest.cmake
  EXCLUDE_FROM_ALL
  )

package_get_option_name(gtest _opt_name)
if(AKANTU_TESTS)
  set(${_opt_name} ON CACHE BOOL "Add GTest support for tests (forced)" FORCE)
endif()

package_add_third_party_script_variable(google-test
  GTEST_VERSION "master")
package_add_third_party_script_variable(google-test
  GTEST_GIT "https://github.com/google/googletest.git")
