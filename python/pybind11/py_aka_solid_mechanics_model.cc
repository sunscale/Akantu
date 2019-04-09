/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
#include "non_linear_solver.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
namespace _aka = akantu;

namespace {

/* -------------------------------------------------------------------------- */

py::module & register_solid_mechanics_models(py::module & mod) {

  py::class_<_aka::SolidMechanicsModelOptions>(mod,
                                               "SolidMechanicsModelOptions")
      .def(py::init<_aka::AnalysisMethod>(),
           py::arg("analysis_method") = _aka::_explicit_lumped_mass);

  py::class_<_aka::Model>(mod, "Model")
      .def("setBaseName", &_aka::Model::setBaseName)
      .def("addDumpFieldVector", &_aka::Model::addDumpFieldVector)
      .def("addDumpField", &_aka::Model::addDumpField)
      .def("dump", &_aka::Model::dump);

  py::class_<_aka::NonLinearSolver>(mod, "NonLinearSolver")
      .def("set",
           [](_aka::NonLinearSolver & self, const std::string & id,
              const _aka::Real & val) {
             if (id == "max_iterations")
               self.set(id, int(val));
             else if (id == "convergence_type")
               self.set(id, akantu::SolveConvergenceCriteria(_aka::UInt(val)));
             else
               self.set(id, val);
           })
      .def("set", &_aka::NonLinearSolver::set<_aka::SolveConvergenceCriteria>);

  py::class_<_aka::ModelSolver>(mod, "ModelSolver")
      .def("getNonLinearSolver",
           (_aka::NonLinearSolver & (_aka::ModelSolver::*)(const _aka::ID &)) &
               _aka::ModelSolver::getNonLinearSolver,
           py::arg("solver_id") = "", py::return_value_policy::reference)
      .def("solveStep", &_aka::ModelSolver::solveStep,
           py::arg("solver_id") = "");

  py::class_<_aka::SolidMechanicsModel, _aka::Model, _aka::ModelSolver>(
      mod, "SolidMechanicsModel")
      .def(py::init<_aka::Mesh &, _aka::UInt, const _aka::ID &,
                    const _aka::MemoryID &, const _aka::ModelType>(),
           py::arg("mesh"),
           py::arg("spatial_dimension") = _aka::_all_dimensions,
           py::arg("id") = "solid_mechanics_model", py::arg("memory_id") = 0,
           py::arg("model_type") = _aka::ModelType::_solid_mechanics_model)
      .def("initFull",
           [](_aka::SolidMechanicsModel & self,
              const _aka::SolidMechanicsModelOptions & options) {
             self.initFull(options);
           },
           py::arg("options") = _aka::SolidMechanicsModelOptions())
      .def("applyBC",
           [](_aka::SolidMechanicsModel & self,
              _aka::BC::DirichletFunctor & func,
              const std::string & element_group) {
             self.applyBC(func, element_group);
           })
      .def("applyBC",
           [](_aka::SolidMechanicsModel & self, _aka::BC::NeumannFunctor & func,
              const std::string & element_group) {
             self.applyBC(func, element_group);
           })
      .def("getExternalForce", &_aka::SolidMechanicsModel::getExternalForce);

  return mod;
} // namespace

} // namespace
