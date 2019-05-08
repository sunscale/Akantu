/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "element_group.hh"
#include "mesh.hh"
#include "py_aka_array.hh"
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */

namespace py = pybind11;

namespace akantu {

/* -------------------------------------------------------------------------- */
__attribute__((visibility("default"))) void register_mesh(py::module & mod) {

  py::class_<ElementGroup>(mod, "ElementGroup");

  py::class_<GroupManager>(mod, "GroupManager")
      .def("getElementGroup",
           [](GroupManager & self, const std::string & name) -> decltype(auto) {
             return self.getElementGroup(name);
           },
           py::return_value_policy::reference);

  py::class_<Mesh, GroupManager>(mod, "Mesh")
      .def(py::init<UInt, const ID &, const MemoryID &>(),
           py::arg("spatial_dimension"), py::arg("id") = "mesh",
           py::arg("memory_id") = 0)
      .def("read", &Mesh::read, py::arg("filename"),
           py::arg("mesh_io_type") = _miot_auto, "read the mesh from a file")
      .def("getNodes",
           [](Mesh & self) -> Array<Real> { return self.getNodes(); },
           py::return_value_policy::reference)
      .def("getNbNodes", &Mesh::getNbNodes)
      .def("distribute", [](Mesh & self) { self.distribute(); });
}
} // namespace akantu
