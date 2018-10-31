#include <pybind11/pybind11.h>
#include <aka_array.hh>

namespace py = pybind11;
using namespace akantu;

PYBIND11_MODULE(py11_akantu_test_common, mod) {
  mod.doc() = "Akantu Test function for common ";

  Array<Real> data{1000, 3};

  mod.def("getData", [&data]() -> Array<Real>& { return data; });
} // Module akantu_test_common
