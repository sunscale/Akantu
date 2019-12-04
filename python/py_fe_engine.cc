/* -------------------------------------------------------------------------- */
#include "py_aka_array.hh"
#include "py_aka_common.hh"
/* -------------------------------------------------------------------------- */
#include <fe_engine.hh>
#include <integration_point.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {

__attribute__((visibility("default"))) void
register_fe_engine(py::module & mod) {

  py::class_<Element>(mod, "Element");

  py::class_<FEEngine>(mod, "FEEngine")
      .def(
          "gradientOnIntegrationPoints",
          [](FEEngine & fem, const Array<Real> & u, Array<Real> & nablauq,
             const UInt nb_degree_of_freedom, const ElementType & type,
             const GhostType & ghost_type,
             const Array<UInt> & filter_elements) {
            fem.gradientOnIntegrationPoints(u, nablauq, nb_degree_of_freedom,
                                            type, ghost_type, filter_elements);
          },
          py::arg("u"), py::arg("nablauq"), py::arg("nb_degree_of_freedom"),
          py::arg("type"), py::arg("ghost_type") = _not_ghost,
          py::arg("filter_elements") = empty_filter)
      .def(
          "interpolateOnIntegrationPoints",
          [](FEEngine & self, const Array<Real> & u, Array<Real> & uq,
             UInt nb_degree_of_freedom, const ElementType & type,
             const GhostType & ghost_type,
             const Array<UInt> & filter_elements) {
            self.interpolateOnIntegrationPoints(
                u, uq, nb_degree_of_freedom, type, ghost_type, filter_elements);
          },
          py::arg("u"), py::arg("uq"), py::arg("nb_degree_of_freedom"),
          py::arg("type"), py::arg("ghost_type") = _not_ghost,
          py::arg("filter_elements") = empty_filter)
      .def(
          "interpolateOnIntegrationPoints",
          [](FEEngine & self, const Array<Real> & u,
             ElementTypeMapArray<Real> & uq,
             const ElementTypeMapArray<UInt> * filter_elements) {
            self.interpolateOnIntegrationPoints(u, uq, filter_elements);
          },
          py::arg("u"), py::arg("uq"), py::arg("filter_elements") = nullptr)
      .def(
          "computeIntegrationPointsCoordinates",
          [](FEEngine & self, ElementTypeMapArray<Real> & coordinates,
             const ElementTypeMapArray<UInt> * filter_elements)
              -> decltype(auto) {
            return self.computeIntegrationPointsCoordinates(coordinates,
                                                            filter_elements);
          },
          py::arg("coordinates"), py::arg("filter_elements") = nullptr)
      .def(
          "assembleFieldLumped",
          [](FEEngine & fem,
             const std::function<void(Matrix<Real> &, const Element &)> &
                 field_funct,
             const ID & matrix_id, const ID & dof_id, DOFManager & dof_manager,
             ElementType type, const GhostType & ghost_type) {
            fem.assembleFieldLumped(field_funct, matrix_id, dof_id, dof_manager,
                                    type, ghost_type);
          },
          py::arg("field_funct"), py::arg("matrix_id"), py::arg("dof_id"),
          py::arg("dof_manager"), py::arg("type"),
          py::arg("ghost_type") = _not_ghost)
      .def(
          "assembleFieldMatrix",
          [](FEEngine & fem,
             const std::function<void(Matrix<Real> &, const Element &)> &
                 field_funct,
             const ID & matrix_id, const ID & dof_id, DOFManager & dof_manager,
             ElementType type, const GhostType & ghost_type = _not_ghost) {
            fem.assembleFieldMatrix(field_funct, matrix_id, dof_id, dof_manager,
                                    type, ghost_type);
          },
          py::arg("field_funct"), py::arg("matrix_id"), py::arg("dof_id"),
          py::arg("dof_manager"), py::arg("type"),
          py::arg("ghost_type") = _not_ghost);

  py::class_<IntegrationPoint>(mod, "IntegrationPoint");
}
} // namespace akantu
