/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "py_aka_array.cc"
/* -------------------------------------------------------------------------- */

namespace py = pybind11;

namespace akantu {

/* -------------------------------------------------------------------------- */
__attribute__((visibility("default"))) void register_mesh(py::module & mod) {
  py::class_<Mesh>(mod, "Mesh")
      .def(py::init<UInt, const ID &, const MemoryID &>(),
           py::arg("spatial_dimension"), py::arg("id") = "mesh",
           py::arg("memory_id") = 0)
      .def("read", &Mesh::read, py::arg("filename"),
           py::arg("mesh_io_type") = _miot_auto, "read the mesh from a file")
      .def("getNodes",
           [](Mesh & self) -> Array<Real> { return self.getNodes(); },
           py::return_value_policy::reference)
      .def("getNbNodes", &Mesh::getNbNodes);
}
} // namespace akantu
