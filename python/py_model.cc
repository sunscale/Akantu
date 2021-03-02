/* -------------------------------------------------------------------------- */
#include "py_aka_array.hh"
/* -------------------------------------------------------------------------- */
#include <model.hh>
#include <non_linear_solver.hh>
#include <sparse_matrix_aij.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
void register_model(py::module & mod) {
  py::class_<DOFManager>(mod, "DOFManager")
      .def("getMatrix", &DOFManager::getMatrix,
           py::return_value_policy::reference)
      .def(
          "getNewMatrix",
          [](DOFManager & self, const std::string & name,
             const std::string & matrix_to_copy_id) -> decltype(auto) {
            return self.getNewMatrix(name, matrix_to_copy_id);
          },
          py::return_value_policy::reference)
      .def(
          "getResidual",
          [](DOFManager & self) -> decltype(auto) {
            return self.getResidual();
          },
          py::return_value_policy::reference)
      .def("getArrayPerDOFs", &DOFManager::getArrayPerDOFs)
      .def("assembleToResidual", &DOFManager::assembleToResidual);

  py::class_<NonLinearSolver>(mod, "NonLinearSolver")
      .def(
          "set",
          [](NonLinearSolver & self, const std::string & id, const Real & val) {
            if (id == "max_iterations") {
              self.set(id, int(val));
            } else {
              self.set(id, val);
            }
          })
      .def("set",
           [](NonLinearSolver & self, const std::string & id,
              const SolveConvergenceCriteria & val) { self.set(id, val); });

  py::class_<ModelSolver, Parsable>(mod, "ModelSolver",
                                    py::multiple_inheritance())
      .def("getNonLinearSolver",
           (NonLinearSolver & (ModelSolver::*)(const ID &)) &
               ModelSolver::getNonLinearSolver,
           py::arg("solver_id") = "", py::return_value_policy::reference)
      .def("solveStep", [](ModelSolver & self) { self.solveStep(); })
      .def("solveStep", [](ModelSolver & self, const ID & solver_id) {
        self.solveStep(solver_id);
      });

  py::class_<Model, ModelSolver>(mod, "Model", py::multiple_inheritance())
      .def("setBaseName", &Model::setBaseName)
      .def("setDirectory", &Model::setDirectory)
      .def("getFEEngine", &Model::getFEEngine, py::arg("name") = "",
           py::return_value_policy::reference)
      .def("getFEEngineBoundary", &Model::getFEEngine, py::arg("name") = "",
           py::return_value_policy::reference)
      .def("addDumpFieldVector", &Model::addDumpFieldVector)
      .def("addDumpField", &Model::addDumpField)
      .def("setBaseNameToDumper", &Model::setBaseNameToDumper)
      .def("addDumpFieldVectorToDumper", &Model::addDumpFieldVectorToDumper)
      .def("addDumpFieldToDumper", &Model::addDumpFieldToDumper)
      .def("dump", &Model::dump)
      .def("initNewSolver", &Model::initNewSolver)
      .def("getNewSolver", [](Model & self, const std::string id, const TimeStepSolverType & time,
			      const NonLinearSolverType & type) {
	     self.getNewSolver(id, time, type);
	   }, py::return_value_policy::reference)
      .def("setIntegrationScheme", [](Model & self, const std::string id, const std::string primal,
				      const IntegrationSchemeType & scheme) {
	     self.setIntegrationScheme(id, primal, scheme);
	   })
      .def("getDOFManager", &Model::getDOFManager,
           py::return_value_policy::reference);
}

} // namespace akantu
