/* -------------------------------------------------------------------------- */
#include "py_aka_error.hh"
/* -------------------------------------------------------------------------- */
#include <aka_error.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {
/* -------------------------------------------------------------------------- */

[[gnu::visibility("default")]] void register_error(py::module & mod) {

  mod.def("setDebugLevel", &debug::setDebugLevel);
  mod.def("getDebugLevel", &debug::getDebugLevel);
  mod.def("printBacktrace", [](bool flag) { debug::debugger.printBacktrace(flag); });

  py::enum_<DebugLevel>(mod, "DebugLevel")
      .value("dblError", dblError)
      .value("dblException", dblException)
      .value("dblCritical", dblCritical)
      .value("dblMajor", dblMajor)
      .value("dblWarning", dblWarning)
      .value("dblInfo", dblInfo)
      .value("dblTrace", dblTrace)
      .value("dblAccessory", dblAccessory)
      .value("dblDebug", dblDebug)
      .value("dblDump", dblDump)
      .value("dblTest", dblTest)
      .export_values();
}

} // namespace akantu
