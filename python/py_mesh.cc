/* -------------------------------------------------------------------------- */
#include "aka_config.hh"
/* -------------------------------------------------------------------------- */
#include "py_aka_array.hh"
/* -------------------------------------------------------------------------- */
#include <mesh.hh>
#include <mesh_utils.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {


/* -------------------------------------------------------------------------- */
template <typename T>
void register_element_type_map_array(py::module & mod,
                                     const std::string & name) {
  py::class_<ElementTypeMapArray<T>, std::shared_ptr<ElementTypeMapArray<T>>>(
      mod, ("ElementTypeMapArray" + name).c_str())
      .def(
          "__call__",
          [](ElementTypeMapArray<T> & self, ElementType type,
             GhostType ghost_type) -> decltype(auto) {
            return self(type, ghost_type);
          },
          py::arg("type"), py::arg("ghost_type") = _not_ghost,
          py::return_value_policy::reference)
      .def(
          "elementTypes",
          [](ElementTypeMapArray<T> & self, UInt _dim, GhostType _ghost_type,
             ElementKind _kind) -> decltype(auto) {
            auto types = self.elementTypes(_dim, _ghost_type, _kind);
            std::vector<ElementType> _types;
            for (auto && t : types) {
              _types.push_back(t);
            }
            return _types;
          },
          py::arg("dim") = _all_dimensions, py::arg("ghost_type") = _not_ghost,
          py::arg("kind") = _ek_regular);
}


/* -------------------------------------------------------------------------- */
void register_mesh(py::module & mod) {

  register_element_type_map_array<Real>(mod, "Real");
  register_element_type_map_array<UInt>(mod, "UInt");

  py::class_<MeshData>(mod, "MeshData")
      .def(
          "getElementalDataUInt",
          [](MeshData & _this, const ID & name) -> ElementTypeMapArray<UInt> & {
            return _this.getElementalData<UInt>(name);
          },
          py::return_value_policy::reference);

  py::class_<Mesh, GroupManager, Dumpable, MeshData>(mod, "Mesh",
                                                     py::multiple_inheritance())
      .def(py::init<UInt, const ID &, const MemoryID &>(),
           py::arg("spatial_dimension"), py::arg("id") = "mesh",
           py::arg("memory_id") = 0)
      .def("read", &Mesh::read, py::arg("filename"),
           py::arg("mesh_io_type") = _miot_auto, "read the mesh from a file")
      .def(
          "getNodes",
          [](Mesh & self) -> decltype(auto) { return self.getNodes(); },
          py::return_value_policy::reference)
      .def("getNbNodes", &Mesh::getNbNodes)
      .def(
          "getConnectivity",
          [](Mesh & self, ElementType type) -> decltype(auto) {
            return self.getConnectivity(type);
          },
          py::return_value_policy::reference)
      .def("distribute", [](Mesh & self) { self.distribute(); })
      .def("makePeriodic",
           [](Mesh & self, const SpatialDirection & direction) {
             self.makePeriodic(direction);
           })
      .def(
          "getNbElement",
          [](Mesh & self, const UInt spatial_dimension,
             GhostType ghost_type, ElementKind kind) {
            return self.getNbElement(spatial_dimension, ghost_type, kind);
          },
          py::arg("spatial_dimension") = _all_dimensions,
          py::arg("ghost_type") = _not_ghost, py::arg("kind") = _ek_not_defined)
      .def(
          "getNbElement",
          [](Mesh & self, ElementType type,
             GhostType ghost_type) {
            return self.getNbElement(type, ghost_type);
          },
          py::arg("type"), py::arg("ghost_type") = _not_ghost)
      .def_static("getSpatialDimension", [](ElementType & type) {
        return Mesh::getSpatialDimension(type);
      });

  /* ------------------------------------------------------------------------ */
  py::class_<MeshUtils>(mod, "MeshUtils")
      .def_static("buildFacets", &MeshUtils::buildFacets);
}
} // namespace akantu
