/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "parameter_registry.hh"
#include "parsable.hh"
#include "parser.hh"
#include <map>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {
std::map<void *, std::map<std::string, void *>> map_params;

__attribute__((visibility("default"))) void register_parser(py::module & mod) {

  py::class_<ParameterRegistry>(mod, "ParameterRegistry",
                                py::multiple_inheritance())
      .def("registerParamReal",
           [](ParameterRegistry & self, const std::string & name, UInt type,
              const std::string & description) {
             Real * p = new Real;
             map_params[&self][name] = p;
             self.registerParam<Real>(name, *p, ParameterAccessType(type),
                                      description);
           })
      .def("registerParamReal",
           [](ParameterRegistry & self, const Real & _default,
              const std::string & name, UInt type,
              const std::string & description) {
             Real * p = new Real;
             map_params[&self][name] = p;
             self.registerParam<Real>(name, *p, _default,
                                      ParameterAccessType(type), description);
           })
      .def("getReal", [](ParameterRegistry & self, const std::string & name) {
        return Real(self.get(name));
      });

  py::class_<Parsable, ParameterRegistry>(mod, "Parsable",
                                          py::multiple_inheritance())
      .def(py::init<const ParserType &, const ID &>());

  mod.def("parseInput",
          [](const std::string & input_file) {
            getStaticParser().parse(input_file);
          },
          "Parse an Akantu input file");
}
} // namespace akantu
