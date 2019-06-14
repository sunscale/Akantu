/* -------------------------------------------------------------------------- */
#include "py_aka_array.hh"
/* -------------------------------------------------------------------------- */
#include <aka_common.hh>
#include <dumper_iohelper_paraview.hh>
#include <element_group.hh>
#include <mesh.hh>
#include <mesh_utils.hh>
/* -------------------------------------------------------------------------- */
#include <dumpable_inline_impl.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */

namespace py = pybind11;

namespace akantu {

/* -------------------------------------------------------------------------- */
__attribute__((visibility("default"))) void register_mesh(py::module & mod) {

  py::module dumper_module("dumper");
  mod.attr("dumper") = dumper_module;

  py::class_<dumper::Field, std::shared_ptr<dumper::Field>>(dumper_module,
                                                            "Field");
  py::class_<dumper::ElementalField<UInt>, dumper::Field,
             std::shared_ptr<dumper::ElementalField<UInt>>>(
      dumper_module, "ElementalFieldUInt", py::multiple_inheritance())
      .def(py::init<dumper::ElementalField<UInt>::field_type &, UInt, GhostType,
                    ElementKind>(),
           py::arg("field"), py::arg("spatial_dimension") = _all_dimensions,
           py::arg("ghost_type") = _not_ghost,
           py::arg("element_kind") = _ek_not_defined);

  py::class_<NodeGroup>(mod, "NodeGroup");

  py::class_<ElementGroup>(mod, "ElementGroup")
      .def("getNodes", &ElementGroup::getNodes,
           py::return_value_policy::reference)
      .def("getNodeGroup",
           [](ElementGroup & self) -> decltype(auto) {
             return self.getNodeGroup();
           },
           py::return_value_policy::reference);

  py::class_<MeshData>(mod, "MeshData")
      .def(
          "getElementalDataUInt",
          [](MeshData & _this, const ID & name) -> ElementTypeMapArray<UInt> & {
            return _this.getElementalData<UInt>(name);
          },
          py::return_value_policy::reference);

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
      .def("addDumpFieldExternal",
           [](Dumpable & _this, const std::string & field_id,
              std::shared_ptr<dumper::Field> field) {
             return _this.addDumpFieldExternal(field_id, field);
           },
           py::arg("field_id"), py::arg("field"))
      .def("addDumpFieldExternalToDumper",
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

  py::class_<GroupManager>(mod, "GroupManager")
      .def("getElementGroup",
           [](GroupManager & self, const std::string & name) -> decltype(auto) {
             return self.getElementGroup(name);
           },
           py::return_value_policy::reference)
      .def("createNodeGroup", &GroupManager::createNodeGroup,
           py::return_value_policy::reference)
      .def("createElementGroup",
           py::overload_cast<const std::string &, UInt, bool>(
               &GroupManager::createElementGroup),
           py::return_value_policy::reference)
      .def("createGroupsFromMeshDataUInt",
           &GroupManager::createGroupsFromMeshData<UInt>)
      .def("createElementGroupFromNodeGroup",
           &GroupManager::createElementGroupFromNodeGroup, py::arg("name"),
           py::arg("node_group"), py::arg("dimension") = _all_dimensions)
      .def("getNodeGroup",
           [](GroupManager & self, const std::string & name) -> decltype(auto) {
             return self.getNodeGroup(name);
           },
           py::return_value_policy::reference)
      .def("createBoundaryGroupFromGeometry",
           &GroupManager::createBoundaryGroupFromGeometry)
      .def("getElementGroups", &GroupManager::getElementGroups,
           py::return_value_policy::reference);

  py::class_<Mesh, GroupManager, Dumpable, MeshData>(mod, "Mesh",
                                                     py::multiple_inheritance())
      .def(py::init<UInt, const ID &, const MemoryID &>(),
           py::arg("spatial_dimension"), py::arg("id") = "mesh",
           py::arg("memory_id") = 0)
      .def("read", &Mesh::read, py::arg("filename"),
           py::arg("mesh_io_type") = _miot_auto, "read the mesh from a file")
      .def("getNodes",
           [](Mesh & self) -> decltype(auto) { return self.getNodes(); },
           py::return_value_policy::reference)
      .def("getNbNodes", &Mesh::getNbNodes)
      .def("distribute", [](Mesh & self) { self.distribute(); })
      .def("getNbElement",
           [](Mesh & self, const UInt spatial_dimension,
              const GhostType & ghost_type, const ElementKind & kind) {
             return self.getNbElement(spatial_dimension, ghost_type, kind);
           },
           py::arg("spatial_dimension") = _all_dimensions,
           py::arg("ghost_type") = _not_ghost,
           py::arg("kind") = _ek_not_defined)
      .def("getNbElement",
           [](Mesh & self, const ElementType & type,
              const GhostType & ghost_type) {
             return self.getNbElement(type, ghost_type);
           },
           py::arg("type"), py::arg("ghost_type") = _not_ghost)
      .def_static("getSpatialDimension", [](ElementType & type) {
        return Mesh::getSpatialDimension(type);
      });

  py::class_<MeshUtils>(mod, "MeshUtils")
      .def_static("buildFacets", &MeshUtils::buildFacets);
}
} // namespace akantu
