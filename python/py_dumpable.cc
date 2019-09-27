/* -------------------------------------------------------------------------- */
#include "py_aka_array.hh"
/* -------------------------------------------------------------------------- */
#include <dumper_iohelper_paraview.hh>
#include <mesh.hh>
/* -------------------------------------------------------------------------- */
#include <dumpable_inline_impl.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {

std::vector<detail::ArrayProxy<Real>> tmp_array;

void register_dumpable(py::module & mod) {
  /* ------------------------------------------------------------------------ */
  py::class_<Dumpable>(mod, "Dumpable")
      .def("registerDumperParaview", &Dumpable::registerDumper<DumperParaview>,
           py::arg("dumper_name"), py::arg("file_name"),
           py::arg("is_default") = false)
      .def("addDumpMeshToDumper", &Dumpable::addDumpMeshToDumper,
           py::arg("dumper_name"), py::arg("mesh"), py::arg("dimension"),
           py::arg("ghost_type") = _not_ghost,
           py::arg("element_kind") = _ek_regular)
      .def("addDumpMesh", &Dumpable::addDumpMesh, py::arg("mesh"),
           py::arg("dimension"), py::arg("ghost_type") = _not_ghost,
           py::arg("element_kind") = _ek_regular)
      .def("addDumpField", &Dumpable::addDumpField, py::arg("field_id"))
      .def("addDumpFieldToDumper", &Dumpable::addDumpFieldToDumper,
           py::arg("dumper_name"), py::arg("field_id"))
      .def(
          "addDumpFieldExternal",
          [](Dumpable & _this, const std::string & field_id,
             std::shared_ptr<dumper::Field> field) {
            return _this.addDumpFieldExternal(field_id, field);
          },
          py::arg("field_id"), py::arg("field"))
      .def(
          "addDumpFieldExternal",
          [](Dumpable & _this, const std::string & field_id,
             Array<Real> & field) {
            auto & tmp = dynamic_cast<detail::ArrayProxy<Real>&>(field);
            tmp_array.push_back(tmp);
            std::cout << tmp.storage() << std::endl;
            return _this.addDumpFieldExternal(field_id, tmp_array.back());
          },
          py::arg("field_id"), py::arg("field"))
      .def(
          "addDumpFieldExternalToDumper",
          [](Dumpable & _this, const std::string & dumper_name,
             const std::string & field_id,
             std::shared_ptr<dumper::Field> field) {
            return _this.addDumpFieldExternalToDumper(dumper_name, field_id,
                                                      field);
          },
          py::arg("dumper_name"), py::arg("field_id"), py::arg("field"))

      .def("dump", py::overload_cast<>(&Dumpable::dump))
      .def("dump", py::overload_cast<Real, UInt>(&Dumpable::dump),
           py::arg("time"), py::arg("step"))
      .def("dump", py::overload_cast<UInt>(&Dumpable::dump), py::arg("step"))
      .def("dump",
           py::overload_cast<const std::string &, UInt>(&Dumpable::dump),
           py::arg("dumper_name"), py::arg("step"))
      .def("dump",
           py::overload_cast<const std::string &, Real, UInt>(&Dumpable::dump),
           py::arg("dumper_name"), py::arg("time"), py::arg("step"))
      .def("dump", py::overload_cast<const std::string &>(&Dumpable::dump),
           py::arg("dumper_name"));

  /* ------------------------------------------------------------------------ */
  py::module dumper_module("dumper");
  mod.attr("dumper") = dumper_module;

  /* ------------------------------------------------------------------------ */
  py::class_<dumper::Field, std::shared_ptr<dumper::Field>>(dumper_module,
                                                            "Field");

  /* ------------------------------------------------------------------------ */
  py::class_<dumper::ElementalField<UInt>, dumper::Field,
             std::shared_ptr<dumper::ElementalField<UInt>>>(
      dumper_module, "ElementalFieldUInt", py::multiple_inheritance())
      .def(py::init<dumper::ElementalField<UInt>::field_type &, UInt, GhostType,
                    ElementKind>(),
           py::arg("field"), py::arg("spatial_dimension") = _all_dimensions,
           py::arg("ghost_type") = _not_ghost,
           py::arg("element_kind") = _ek_not_defined);
}

} // namespace akantu
