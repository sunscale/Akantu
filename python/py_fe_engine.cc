/* -------------------------------------------------------------------------- */
#include <fe_engine.hh>
#include <integration_point.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {

__attribute__((visibility("default"))) void
register_fe_engine(py::module & mod) {

  py::class_<Element>(mod, "Element");

  py::class_<FEEngine>(mod, "FEEngine")
      .def("gradientOnIntegrationPoints",
           [](FEEngine & fem, const Array<Real> & u, Array<Real> & nablauq,
              const UInt nb_degree_of_freedom, const ElementType & type,
              const GhostType & ghost_type = _not_ghost,
              const Array<UInt> & filter_elements = empty_filter) {
             fem.gradientOnIntegrationPoints(u, nablauq, nb_degree_of_freedom,
                                             type, ghost_type, filter_elements);
           })
      .def("interpolateOnIntegrationPoints",
           [](FEEngine & self, const Array<Real> & u, Array<Real> & uq,
              UInt nb_degree_of_freedom, const ElementType & type,
              const GhostType & ghost_type = _not_ghost,
              const Array<UInt> & filter_elements = empty_filter) {
             self.interpolateOnIntegrationPoints(u, uq, nb_degree_of_freedom,
                                                 type, ghost_type,
                                                 filter_elements);
           })
      .def("interpolateOnIntegrationPoints",
           [](FEEngine & self, const Array<Real> & u,
              ElementTypeMapArray<Real> & uq,
              const ElementTypeMapArray<UInt> * filter_elements = nullptr) {
             self.interpolateOnIntegrationPoints(u, uq, filter_elements);
           })
      .def("computeIntegrationPointsCoordinates",
           [](FEEngine & self, ElementTypeMapArray<Real> & coordinates,
              const ElementTypeMapArray<UInt> * filter_elements =
                  nullptr) -> decltype(auto) {
             return self.computeIntegrationPointsCoordinates(coordinates,
                                                             filter_elements);
           })
      .def("assembleFieldLumped",
           [](FEEngine & fem,
              const std::function<void(Matrix<Real> &, const Element &)> &
                  field_funct,
              const ID & matrix_id, const ID & dof_id, DOFManager & dof_manager,
              ElementType type, const GhostType & ghost_type = _not_ghost) {
             fem.assembleFieldLumped(field_funct, matrix_id, dof_id,
                                     dof_manager, type, ghost_type);
           })
      .def("assembleFieldMatrix",
           [](FEEngine & fem,
              const std::function<void(Matrix<Real> &, const Element &)> &
                  field_funct,
              const ID & matrix_id, const ID & dof_id, DOFManager & dof_manager,
              ElementType type, const GhostType & ghost_type = _not_ghost) {
             fem.assembleFieldMatrix(field_funct, matrix_id, dof_id,
                                     dof_manager, type, ghost_type);
           });

  py::class_<IntegrationPoint>(mod, "IntegrationPoint");
}
} // namespace akantu
