/* -------------------------------------------------------------------------- */
#include "py_aka_array.hh"
/* -------------------------------------------------------------------------- */
#include <non_linear_solver.hh>
#include <solid_mechanics_model.hh>
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
  def(                                                                         \
      #func_name,                                                              \
      [](SolidMechanicsModel & self) -> decltype(auto) {                       \
        return self.func_name();                                               \
      },                                                                       \
      py::return_value_policy::reference)

#define def_function(func_name)                                                \
  def(#func_name, [](SolidMechanicsModel & self) -> decltype(auto) {           \
    return self.func_name();                                                   \
  })
/* -------------------------------------------------------------------------- */

void
register_solid_mechanics_model(py::module & mod) {

  py::class_<SolidMechanicsModelOptions>(mod, "SolidMechanicsModelOptions")
      .def(py::init<AnalysisMethod>(),
           py::arg("_analysis_method") = _explicit_lumped_mass);

  py::class_<SolidMechanicsModel, Model>(mod, "SolidMechanicsModel",
                                         py::multiple_inheritance())
      .def(py::init<Mesh &, UInt, const ID &,
                    const ModelType>(),
           py::arg("mesh"), py::arg("spatial_dimension") = _all_dimensions,
           py::arg("id") = "solid_mechanics_model",
           py::arg("model_type") = ModelType::_solid_mechanics_model)
      .def(
          "initFull",
          [](SolidMechanicsModel & self,
             const SolidMechanicsModelOptions & options) {
            self.initFull(options);
          },
          py::arg("option") = SolidMechanicsModelOptions())
      .def(
          "initFull",
          [](SolidMechanicsModel & self,
             const AnalysisMethod & analysis_method) {
            self.initFull(_analysis_method = analysis_method);
          },
          py::arg("_analysis_method"))
      .def_deprecated("applyDirichletBC", "Deprecated: use applyBC")
      .def("applyBC",
           [](SolidMechanicsModel & self,
              BC::Dirichlet::DirichletFunctor & func,
              const std::string & element_group) {
             self.applyBC(func, element_group);
           })
      .def("applyBC",
           [](SolidMechanicsModel & self, BC::Neumann::NeumannFunctor & func,
              const std::string & element_group) {
             self.applyBC(func, element_group);
           })
      .def("setTimeStep", &SolidMechanicsModel::setTimeStep,
           py::arg("time_step"), py::arg("solver_id") = "")
      .def("getEnergy",
           py::overload_cast<const std::string &>(
               &SolidMechanicsModel::getEnergy),
           py::arg("energy_id"))
      .def_function(assembleStiffnessMatrix)
      .def_function(assembleInternalForces)
      .def_function(assembleMass)
      .def_function(assembleMassLumped)
      .def_function(getStableTimeStep)
      .def_function_nocopy(getExternalForce)
      .def_function_nocopy(getDisplacement)
      .def_function_nocopy(getPreviousDisplacement)
      .def_function_nocopy(getIncrement)
      .def_function_nocopy(getInternalForce)
      .def_function_nocopy(getMass)
      .def_function_nocopy(getVelocity)
      .def_function_nocopy(getAcceleration)
      .def_function_nocopy(getInternalForce)
      .def_function_nocopy(getBlockedDOFs)
      .def_function_nocopy(getMesh)
      .def("dump", py::overload_cast<>(&SolidMechanicsModel::dump))
      .def("dump",
           py::overload_cast<const std::string &>(&SolidMechanicsModel::dump))
      .def("dump", py::overload_cast<const std::string &, UInt>(
                       &SolidMechanicsModel::dump))
      .def("dump", py::overload_cast<const std::string &, Real, UInt>(
                       &SolidMechanicsModel::dump))
      .def("getMaterial",
           py::overload_cast<UInt>(&SolidMechanicsModel::getMaterial),
           py::return_value_policy::reference)
      .def("getMaterial",
           py::overload_cast<const std::string &>(
               &SolidMechanicsModel::getMaterial),
           py::return_value_policy::reference)
      .def("getMaterialIndex", &SolidMechanicsModel::getMaterialIndex)
      .def("setMaterialSelector", &SolidMechanicsModel::setMaterialSelector)
      .def("getMaterialSelector", &SolidMechanicsModel::getMaterialSelector);
}

} // namespace akantu
