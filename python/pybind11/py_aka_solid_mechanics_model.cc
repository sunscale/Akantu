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

/* -------------------------------------------------------------------------- */

namespace akantu {

class MaterialPython : public Material {

public:
  MaterialPython(SolidMechanicsModel & model, PyObject * obj,
                 const ID & id = "");

  ~MaterialPython() override = default;

public:
  void registerInternals();

  void initMaterial() override;

  void computeStress(ElementType el_type,
                     GhostType ghost_type = _not_ghost) override;

  void computeTangentModuli(const ElementType & el_type,
                            Array<Real> & tangent_matrix,
                            GhostType ghost_type = _not_ghost) override;

  Real getPushWaveSpeed(const Element & element) const override;

  Real getEnergy(const std::string & type) override;

  Real getEnergyForType(const std::string & type, ElementType el_type);

protected:
  std::map<std::string, Real> local_params;
  std::map<std::string, std::unique_ptr<InternalField<Real>>> internals;
};

} // namespace akantu

/* -------------------------------------------------------------------------- */
namespace _aka = akantu;

namespace {

/* -------------------------------------------------------------------------- */
#define def_deprecated(func_name, mesg)                                        \
  def(func_name, [](py::args, py::kwargs) { AKANTU_ERROR(mesg); })

#define def_function_nocopy(func_name)                                         \
  def(#func_name,                                                              \
      [](_aka::SolidMechanicsModel & self) -> decltype(auto) {                 \
        return self.func_name();                                               \
      },                                                                       \
      py::return_value_policy::reference)

#define def_function(func_name)                                                \
  def(#func_name, [](_aka::SolidMechanicsModel & self) -> decltype(auto) {     \
    return self.func_name();                                                   \
  })
/* -------------------------------------------------------------------------- */

py::module & register_solid_mechanics_models(py::module & mod) {

  py::class_<_aka::SolidMechanicsModelOptions>(
      mod, "_aka::SolidMechanicsModelOptions")
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
      mod, "_aka::SolidMechanicsModel")
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
           py::arg("_analysis_method") = _aka::SolidMechanicsModelOptions())
      .def("initFull",
           [](_aka::SolidMechanicsModel & self,
              const _aka::AnalysisMethod & _analysis_method) {
             self.initFull(_aka::SolidMechanicsModelOptions(_analysis_method));
           },
           py::arg("_analysis_method"))
      .def_deprecated("applyDirichletBC", "Deprecated: use applyBC")
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
      .def("setTimeStep", &_aka::SolidMechanicsModel::setTimeStep,
           py::arg("time_step"), py::arg("solver_id") = "")
      .def("getEnergy",
           py::overload_cast<const std::string &>(
               &_aka::SolidMechanicsModel::getEnergy),
           py::arg("energy_id"))
      .def_function(assembleStiffnessMatrix)
      .def_function(assembleInternalForces)
      .def_function(getStableTimeStep)
      .def_function_nocopy(getExternalForce)
      .def_function_nocopy(getDisplacement)
      .def_function_nocopy(getPreviousDisplacement)
      .def_function_nocopy(getIncrement)
      .def_function_nocopy(getMass)
      .def_function_nocopy(getVelocity)
      .def_function_nocopy(getAcceleration)
      .def_function_nocopy(getInternalForce)
      .def_function_nocopy(getBlockedDOFs)
      .def_function_nocopy(getIncrementFlag);

  py::class_<_aka::MaterialFactory>(mod, "MaterialFactory");

  return mod;
} // namespace

} // namespace
