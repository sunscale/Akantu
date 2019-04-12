#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

#include "py_aka_common.cc"
#include "py_aka_array.cc"
#include "py_aka_boundary_conditions.cc"
#include "py_aka_solid_mechanics_model.cc"

PYBIND11_MODULE(py11_akantu, mod) {
  mod.doc() = "Akantu python interface";

  register_enums(mod);
  register_boundary_conditions(mod);
  register_solid_mechanics_models(mod);

} // Module akantu
