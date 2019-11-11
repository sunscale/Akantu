set(PYBIND11_PYTHON_VERSION ${AKANTU_PREFERRED_PYTHON_VERSION} CACHE INTERNAL "")

package_declare(pybind11 EXTERNAL
  EXTRA_PACKAGE_OPTIONS ARGS "2.4.2;CONFIG" LINK_LIBRARIES pybind11::embed PREFIX pybind11
  DESCRIPTION "Akantu's pybind11 interface"
  SYSTEM AUTO third-party/cmake/pybind11.cmake
  EXCLUDE_FROM_ALL
  )

package_add_third_party_script_variable(pybind11
  PYBIND11_VERSION "v2.4.2")
package_add_third_party_script_variable(pybind11
  PYBIND11_GIT "https://github.com/pybind/pybind11.git")
