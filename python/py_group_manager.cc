/* -------------------------------------------------------------------------- */
#include "py_aka_array.hh"
/* -------------------------------------------------------------------------- */
#include <element_group.hh>
#include <node_group.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
void register_group_manager(py::module & mod) {
  /* ------------------------------------------------------------------------ */
  py::class_<NodeGroup>(mod, "NodeGroup")
      .def(
          "getNodes",
          [](NodeGroup & self) -> decltype(auto) { return self.getNodes(); },
          py::return_value_policy::reference)
      .def("getName", &NodeGroup::getName);

  /* ------------------------------------------------------------------------ */
  py::class_<ElementGroup>(mod, "ElementGroup")
      .def(
          "getNodeGroup",
          [](ElementGroup & self) -> decltype(auto) {
            return self.getNodeGroup();
          },
          py::return_value_policy::reference)
      .def("getName", &ElementGroup::getName)
      .def(
          "getElements",
          [](ElementGroup & self) -> decltype(auto) {
            return self.getElements();
          },
          py::return_value_policy::reference);

  /* ------------------------------------------------------------------------ */
  py::class_<GroupManager>(mod, "GroupManager")
      .def(
          "getElementGroup",
          [](GroupManager & self, const std::string & name) -> decltype(auto) {
            return self.getElementGroup(name);
          },
          py::return_value_policy::reference)
      .def("iterateElementGroups",
           [](GroupManager & self) -> decltype(auto) {
             std::vector<std::reference_wrapper<ElementGroup>> groups;
             for (auto & group : self.iterateElementGroups()) {
               groups.emplace_back(group);
             }
             return groups;
           })
      .def("iterateNodeGroups",
           [](GroupManager & self) -> decltype(auto) {
             std::vector<std::reference_wrapper<NodeGroup>> groups;
             for (auto & group : self.iterateNodeGroups()) {
               groups.emplace_back(group);
             }
             return groups;
           })
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
      .def(
          "getNodeGroup",
          [](GroupManager & self, const std::string & name) -> decltype(auto) {
            return self.getNodeGroup(name);
          },
          py::return_value_policy::reference)
      .def(
          "nodeGroups",
          [](GroupManager & self) {
            std::vector<NodeGroup *> groups;
            for (auto & g : self.iterateNodeGroups()) {
              groups.push_back(&g);
            }
            return groups;
          },
          py::return_value_policy::reference)
      .def(
          "elementGroups",
          [](GroupManager & self) {
            std::vector<ElementGroup *> groups;
            for (auto & g : self.iterateElementGroups()) {
              groups.push_back(&g);
            }
            return groups;
          },
          py::return_value_policy::reference)
      .def("createBoundaryGroupFromGeometry",
           &GroupManager::createBoundaryGroupFromGeometry);
}

} // namespace akantu
