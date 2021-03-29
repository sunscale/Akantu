/* -------------------------------------------------------------------------- */
#include "aka_config.hh"
/* -------------------------------------------------------------------------- */
#include "py_aka_common.hh"
#include "py_aka_error.hh"
#include "py_boundary_conditions.hh"
#include "py_fe_engine.hh"
#include "py_group_manager.hh"
#include "py_mesh.hh"
#include "py_model.hh"
#include "py_parser.hh"
#include "py_solver.hh"

#if defined(AKANTU_USE_IOHELPER)
#include "py_dumpable.hh"
#endif

#if defined(AKANTU_SOLID_MECHANICS)
#include "py_material.hh"
#include "py_solid_mechanics_model.hh"
#endif

#if defined(AKANTU_HEAT_TRANSFER)
#include "py_heat_transfer_model.hh"
#endif

#if defined(AKANTU_COHESIVE_ELEMENT)
#include "py_fragment_manager.hh"
#include "py_solid_mechanics_model_cohesive.hh"
#endif

#if defined(AKANTU_STRUCTURAL_MECHANICS)
#include "py_structural_mechanics_model.hh"
#endif
/* -------------------------------------------------------------------------- */
#include <aka_error.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */

namespace py = pybind11;

namespace akantu {
void register_all(pybind11::module & mod) {
  register_initialize(mod);
  register_enums(mod);
  register_error(mod);
  register_functions(mod);
  register_parser(mod);
  register_solvers(mod);

  register_group_manager(mod);
#if defined(AKANTU_USE_IOHELPER)
  register_dumpable(mod);
#endif
  register_mesh(mod);

  register_fe_engine(mod);

  register_boundary_conditions(mod);
  register_model(mod);
#if defined(AKANTU_HEAT_TRANSFER)
  register_heat_transfer_model(mod);
#endif

#if defined(AKANTU_SOLID_MECHANICS)
  register_solid_mechanics_model(mod);
  register_material(mod);
#endif

#if defined(AKANTU_COHESIVE_ELEMENT)
  register_solid_mechanics_model_cohesive(mod);
  register_fragment_manager(mod);
#endif

#if defined(AKANTU_STRUCTURAL_MECHANICS)
  register_structural_mechanics_model(mod);
#endif
}
} // namespace akantu

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
PYBIND11_MODULE(py11_akantu, mod) {
  mod.doc() = "Akantu python interface";

  static py::exception<akantu::debug::Exception> akantu_exception(mod,
                                                                  "Exception");

  py::register_exception_translator([](std::exception_ptr ptr) {
    try {
      if (ptr) {
        std::rethrow_exception(ptr);
      }
    } catch (akantu::debug::Exception & e) {
      if (akantu::debug::debugger.printBacktrace()) {
        akantu::debug::printBacktrace();
      }
      akantu_exception(e.info().c_str());
    }
  });

  akantu::register_all(mod);

  mod.def("has_mpi", []() {
#if defined(AKANTU_USE_MPI)
    return true;
#else
    return false;
#endif
  });

} // Module akantu
