package_declare(gbenchmark EXTERNAL
  DESCRIPTION "Add Google Benchmark support"
  SYSTEM AUTO third-party/cmake/gbenchmark.cmake
  EXCLUDE_FROM_ALL
  )

package_get_option_name(gbenchmark _opt_name)
package_add_third_party_script_variable(gbenchmark
  GBENCHMARK_VERSION "master")
package_add_third_party_script_variable(gbenchmark
  GBENCHMARK_GIT "https://github.com/google/benchmark.git")
