/* -------------------------------------------------------------------------- */
#include "py_aka_common.hh"
/* -------------------------------------------------------------------------- */
#include <aka_error.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
/* -------------------------------------------------------------------------- */

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

  akantu::register_all(mod);
} // Module akantu
