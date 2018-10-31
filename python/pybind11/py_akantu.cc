#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "py_aka_common.cc"
#include "py_aka_array.cc"

PYBIND11_MODULE(py11_akantu, mod) {
  mod.doc() = "Akantu python interface";

  register_enums(mod);
} // Module akantu
