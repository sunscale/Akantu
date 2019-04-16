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

std::map<void *, std::map<std::string, void *>> map_params;

/* -------------------------------------------------------------------------- */

namespace akantu {

class PyMaterial : public Material {

public:
  /* Inherit the constructors */
  using Material::Material;

  virtual ~PyMaterial(){};
  // void registerInternals();
  // void initMaterial() override;
  // void computeStress(ElementType el_type,
  //                    GhostType ghost_type = _not_ghost) override;
  // void computeTangentModuli(const ElementType & el_type,
  //                           Array<Real> & tangent_matrix,
  //                           GhostType ghost_type = _not_ghost) override;
  // Real getPushWaveSpeed(const Element & element) const override;
  // Real getEnergy(const std::string & type) override;
  // Real getEnergyForType(const std::string & type, ElementType el_type);

  // protected:
  //   std::map<std::string, Real> local_params;
  //   std::map<std::string, std::unique_ptr<InternalField<Real>>> internals;
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

  py::class_<_aka::ParameterRegistry>(mod, "ParameterRegistry")
      .def("registerParamReal",
           [](_aka::ParameterRegistry & self, std::string name,
              _aka::ParameterAccessType type, const std::string & description) {
             _aka::Real * p = new _aka::Real;
             map_params[&self][name] = p;
             self.registerParam<_aka::Real>(name, *p, type, description);
           })
      .def("registerParamReal",
           [](_aka::ParameterRegistry & self, const _aka::Real & _default,
              std::string name, _aka::ParameterAccessType type,
              const std::string & description) {
             _aka::Real * p = new _aka::Real;
             map_params[&self][name] = p;
             self.registerParam<_aka::Real>(name, *p, _default, type,
                                            description);
           });

  py::class_<_aka::Parsable, _aka::ParameterRegistry>(mod, "Parsable")
      .def(py::init<const _aka::ParserType &, const _aka::ID &>());

  py::class_<_aka::Material, _aka::PyMaterial, _aka::Parsable>(mod, "Material")
      .def(py::init<_aka::SolidMechanicsModel &, const _aka::ID &>());

  py::class_<_aka::MaterialFactory>(mod, "MaterialFactory")
      .def_static("getInstance",
                  []() -> _aka::MaterialFactory & {
                    return _aka::Material::getFactory();
                  },
                  py::return_value_policy::reference)
      .def("registerAllocator", [](_aka::MaterialFactory & self,
                                   const std::string id, py::function func) {
        self.registerAllocator(
            id,
            [func, id](_aka::UInt dim, const _aka::ID &,
                       _aka::SolidMechanicsModel & model,
                       const _aka::ID & id) -> std::unique_ptr<_aka::Material> {
              py::object obj = func(dim, id, model, id);
              auto & ptr = py::cast<_aka::Material &>(obj);

              obj.release();
              return std::unique_ptr<_aka::Material>(&ptr);
            });
      });

  return mod;
}

} // namespace
