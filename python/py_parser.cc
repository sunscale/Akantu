/* -------------------------------------------------------------------------- */
#include "py_aka_array.hh"
/* -------------------------------------------------------------------------- */
#include <aka_common.hh>
#include <parameter_registry.hh>
#include <parsable.hh>
#include <parser.hh>
/* -------------------------------------------------------------------------- */
#include <map>
#include <pybind11/pybind11.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {
std::map<void *, std::map<std::string, void *>> map_params;

void register_parser(py::module & mod) {
  py::enum_<ParameterAccessType>(mod, "ParameterAccessType", py::arithmetic())
      .value("_pat_internal", _pat_internal)
      .value("_pat_writable", _pat_writable)
      .value("_pat_readable", _pat_readable)
      .value("_pat_modifiable", _pat_modifiable)
      .value("_pat_parsable", _pat_parsable)
      .value("_pat_parsmod", _pat_parsmod)
      .export_values();

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
      .def("getReal",
           [](ParameterRegistry & self, const std::string & name) {
             return Real(self.get(name));
           })
      .def("getMatrix",
           [](ParameterRegistry & self, const std::string & name) {
             const Matrix<Real> & res =
                 static_cast<const Matrix<Real> &>(self.get(name));
             return res;
           },
           py::return_value_policy::copy);

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
