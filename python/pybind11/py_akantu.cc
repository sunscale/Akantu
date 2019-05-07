/* -------------------------------------------------------------------------- */
#include <aka_error.hh>
/* -------------------------------------------------------------------------- */
#include "py_aka_boundary_conditions.hh"
#include "py_aka_common.hh"
#include "py_aka_error.hh"
#include "py_aka_fe_engine.hh"
#include "py_aka_material.hh"
#include "py_aka_mesh.hh"
#include "py_aka_model.hh"
#include "py_aka_parser.hh"
#include "py_aka_solid_mechanics_model.hh"
#include "py_aka_solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(py11_akantu, mod) {
  mod.doc() = "Akantu python interface";

  static py::exception<akantu::debug::Exception> akantu_exception(mod,
                                                                  "Exception");
  py::register_exception_translator([](std::exception_ptr p) {
    try {
      if (p)
        std::rethrow_exception(p);
    } catch (akantu::debug::Exception & e) {
      if (akantu::debug::debugger.printBacktrace())
        akantu::debug::printBacktrace(15);
      akantu_exception(e.info().c_str());
    }
  });

  akantu::register_initialize(mod);
  akantu::register_enums(mod);
  akantu::register_error(mod);
  akantu::register_parser(mod);
  akantu::register_boundary_conditions(mod);
  akantu::register_fe_engine(mod);
  akantu::register_model(mod);
  akantu::register_solid_mechanics_model(mod);
  akantu::register_solid_mechanics_model_cohesive(mod);
  akantu::register_material(mod);
  akantu::register_mesh(mod);
} // Module akantu
