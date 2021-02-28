/* -------------------------------------------------------------------------- */
#include "py_aka_array.hh"
/* -------------------------------------------------------------------------- */
#include <non_linear_solver.hh>
#include <phase_field_model.hh>
#include <coupler_solid_phasefield.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
#define def_deprecated(func_name, mesg)                                        \
  def(func_name, [](py::args, py::kwargs) { AKANTU_ERROR(mesg); })

#define def_function_nocopy(func_name)                                         \
  def(#func_name,                                                              \
      [](PhaseFieldModel & self) -> decltype(auto) {                           \
        return self.func_name();                                               \
      },                                                                       \
      py::return_value_policy::reference)

#define def_function(func_name)                                                \
  def(#func_name, [](PhaseFieldModel & self) -> decltype(auto) {               \
    return self.func_name();                                                   \
  })
/* -------------------------------------------------------------------------- */


[[gnu::visibility("default")]] void
register_phase_field_model(py::module & mod) {

  py::class_<PhaseFieldModelOptions>(mod, "PhaseFieldModelOptions")
    .def(py::init<AnalysisMethod>(),
	 py::arg("analysis_method") = _static);

  
  py::class_<PhaseFieldModel, Model>(mod, "PhaseFieldModel",
				     py::multiple_inheritance())
      .def(py::init<Mesh &, UInt, const ID &, const MemoryID &,
                    const ModelType>(),
           py::arg("mesh"), py::arg("spatial_dimension") = _all_dimensions,
           py::arg("id") = "phase_field_model", py::arg("memory_id") = 0,
           py::arg("model_type") = ModelType::_phase_field_model)
      .def("initFull",
           [](PhaseFieldModel & self,
              const PhaseFieldModelOptions & options) {
             self.initFull(options);
           },
           py::arg("_analysis_method") = PhaseFieldModelOptions())
      .def("initFull",
           [](PhaseFieldModel & self,
              const AnalysisMethod & analysis_method) {
             self.initFull(_analysis_method = analysis_method);
           },
           py::arg("_analysis_method"))
      .def_deprecated("applyDirichletBC", "Deprecated: use applyBC")
      .def("applyBC",
           [](PhaseFieldModel & self,
              BC::Dirichlet::DirichletFunctor & func,
              const std::string & element_group) {
             self.applyBC(func, element_group);
           })
      .def("applyBC",
           [](PhaseFieldModel & self, BC::Neumann::NeumannFunctor & func,
              const std::string & element_group) {
             self.applyBC(func, element_group);
           })
      .def("setTimeStep", &PhaseFieldModel::setTimeStep,
           py::arg("time_step"), py::arg("solver_id") = "")
      .def_function(assembleStiffnessMatrix)
      .def_function(assembleInternalForces)
      .def_function_nocopy(getDamage)
      .def_function_nocopy(getInternalForce)
      .def_function_nocopy(getBlockedDOFs)
      .def_function_nocopy(getMesh)
      .def("dump", py::overload_cast<>(&PhaseFieldModel::dump))
      .def("dump",
           py::overload_cast<const std::string &>(&PhaseFieldModel::dump))
      .def("dump", py::overload_cast<const std::string &, UInt>(
                       &PhaseFieldModel::dump))
      .def("dump", py::overload_cast<const std::string &, Real, UInt>(
                       &PhaseFieldModel::dump))
      .def("getPhaseField",
           py::overload_cast<UInt>(&PhaseFieldModel::getPhaseField),
           py::return_value_policy::reference)
      .def("getPhaseField",
           py::overload_cast<const std::string &>(
	       &PhaseFieldModel::getPhaseField),
           py::return_value_policy::reference)
      .def("getPhaseFieldIndex", &PhaseFieldModel::getPhaseFieldIndex)
      .def("setPhaseFieldSelector", &PhaseFieldModel::setPhaseFieldSelector);

}
  

[[gnu::visibility("default")]] void
register_phase_field_coupler(py::module & mod) {
  
  py::class_<CouplerSolidPhaseField, Model>(mod, "CouplerSolidPhaseField")
      .def(py::init<Mesh &, UInt, const ID &, const MemoryID &,
	   const ModelType>(),
           py::arg("mesh"), py::arg("spatial_dimension") = _all_dimensions,
           py::arg("id") = "coupler_solid_phasefield", py::arg("memory_id") = 0,
	   py::arg("model_type") = ModelType::_coupler_solid_phasefield)
      .def("solve", [](CouplerSolidPhaseField & self) { self.solve(); })
      .def("getSolidMechanicsModel",
           &CouplerSolidPhaseField::getSolidMechanicsModel,
           py::return_value_policy::reference)
      .def("getPhaseFieldModel",
           &CouplerSolidPhaseField::getPhaseFieldModel,
           py::return_value_policy::reference);

}
  

  
}
