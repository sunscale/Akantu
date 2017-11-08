package_declare(gtest EXTERNAL
  DESCRIPTION "Add GTest support for tests"
  SYSTEM ON third-party/cmake/gtest.cmake
  FALLBACK_ON_SYSTEM_MISSING ON
  )

package_get_option_name(gtest _opt_name)
if(AKANTU_TESTS)
  set(${_opt_name} ON CACHE BOOL "Add GTest support for tests (forced)" FORCE)
endif()

package_add_third_party_script_variable(google-test
  GTEST_VERSION "master")
package_add_third_party_script_variable(google-test
  GTEST_GIT "https://github.com/google/googletest.git")
