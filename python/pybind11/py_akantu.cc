#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "py_aka_boundary_conditions.hh"
#include "py_aka_common.hh"
#include "py_aka_material.hh"
#include "py_aka_mesh.hh"
#include "py_aka_solid_mechanics_model.hh"

PYBIND11_MODULE(py11_akantu, mod) {
  mod.doc() = "Akantu python interface";

  akantu::register_enums(mod);
  akantu::register_boundary_conditions(mod);
  akantu::register_solid_mechanics_model(mod);
  akantu::register_material(mod);
  akantu::register_mesh(mod);

} // Module akantu
