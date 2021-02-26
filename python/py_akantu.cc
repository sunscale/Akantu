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
#include "py_solid_mechanics_model_cohesive.hh"
#endif


#if defined(AKANTU_PHASE_FIELD)
#include "py_phase_field_model.hh"
#endif

/* -------------------------------------------------------------------------- */
#include <aka_error.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
/* -------------------------------------------------------------------------- */

namespace py = pybind11;

namespace akantu {
void register_all(pybind11::module & mod) {
  register_initialize(mod);
  register_enums(mod);
  register_error(mod);
  register_functions(mod);
  register_parser(mod);

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
#endif

#if defined(AKANTU_PHASE_FIELD)
  register_phase_field_model(mod);
#endif
}
} // namespace akantu

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
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

  akantu::register_all(mod);
} // Module akantu
