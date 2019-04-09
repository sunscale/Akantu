/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
#include "mesh.hh"
#include "parser.hh"
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
namespace _aka = akantu;

namespace {

/* -------------------------------------------------------------------------- */

py::module & register_enums(py::module & mod) {
  py::enum_<_aka::SpatialDirection>(mod, "SpatialDirection")
      .value("_x", _aka::_x)
      .value("_y", _aka::_y)
      .value("_z", _aka::_z)
      .export_values();

  py::enum_<_aka::AnalysisMethod>(mod, "AnalysisMethod")
      .value("_static", _aka::_static)
      .value("_implicit_dynamic", _aka::_implicit_dynamic)
      .value("_explicit_lumped_mass", _aka::_explicit_lumped_mass)
      .value("_explicit_lumped_capacity", _aka::_explicit_lumped_capacity)
      .value("_explicit_consistent_mass", _aka::_explicit_consistent_mass)
      .export_values();

  py::enum_<_aka::NonLinearSolverType>(mod, "NonLinearSolverType")
      .value("_nls_linear", _aka::_nls_linear)
      .value("_nls_newton_raphson", _aka::_nls_newton_raphson)
      .value("_nls_newton_raphson_modified", _aka::_nls_newton_raphson_modified)
      .value("_nls_lumped", _aka::_nls_lumped)
      .value("_nls_auto", _aka::_nls_auto)
      .export_values();

  py::enum_<_aka::TimeStepSolverType>(mod, "TimeStepSolverType")
      .value("_tsst_static", _aka::_tsst_static)
      .value("_tsst_dynamic", _aka::_tsst_dynamic)
      .value("_tsst_dynamic_lumped", _aka::_tsst_dynamic_lumped)
      .value("_tsst_not_defined", _aka::_tsst_not_defined)
      .export_values();

  py::enum_<_aka::IntegrationSchemeType>(mod, "IntegrationSchemeType")
      .value("_ist_pseudo_time", _aka::_ist_pseudo_time)
      .value("_ist_forward_euler", _aka::_ist_forward_euler)
      .value("_ist_trapezoidal_rule_1", _aka::_ist_trapezoidal_rule_1)
      .value("_ist_backward_euler", _aka::_ist_backward_euler)
      .value("_ist_central_difference", _aka::_ist_central_difference)
      .value("_ist_fox_goodwin", _aka::_ist_fox_goodwin)
      .value("_ist_trapezoidal_rule_2", _aka::_ist_trapezoidal_rule_2)
      .value("_ist_linear_acceleration", _aka::_ist_linear_acceleration)
      .value("_ist_newmark_beta", _aka::_ist_newmark_beta)
      .value("_ist_generalized_trapezoidal", _aka::_ist_generalized_trapezoidal)
      .export_values();

  py::enum_<_aka::SolveConvergenceCriteria>(mod, "SolveConvergenceCriteria")
      .value("_scc_residual", _aka::_scc_residual)
      .value("_scc_solution", _aka::_scc_solution)
      .value("_scc_residual_mass_wgh", _aka::_scc_residual_mass_wgh)
      .export_values();

  py::enum_<_aka::CohesiveMethod>(mod, "CohesiveMethod")
      .value("_intrinsic", _aka::_intrinsic)
      .value("_extrinsic", _aka::_extrinsic)
      .export_values();

  py::enum_<_aka::GhostType>(mod, "GhostType")
      .value("_not_ghost", _aka::_not_ghost)
      .value("_ghost", _aka::_ghost)
      .value("_casper", _aka::_casper)
      .export_values();

  py::enum_<_aka::MeshIOType>(mod, "MeshIOType")
      .value("_miot_auto", _aka::_miot_auto)
      .value("_miot_gmsh", _aka::_miot_gmsh)
      .value("_miot_gmsh_struct", _aka::_miot_gmsh_struct)
      .value("_miot_diana", _aka::_miot_diana)
      .value("_miot_abaqus", _aka::_miot_abaqus)
      .export_values();

  py::enum_<_aka::ModelType>(mod, "ModelType")

      .value("_model", _aka::ModelType::_model)
      .value("_solid_mechanics_model", _aka::ModelType::_solid_mechanics_model)
      .value("_solid_mechanics_model_cohesive",
             _aka::ModelType::_solid_mechanics_model_cohesive)
      .value("_heat_transfer_model", _aka::ModelType::_heat_transfer_model)
      .value("_structural_mechanics_model",
             _aka::ModelType::_structural_mechanics_model)
      .value("_embedded_model", _aka::ModelType::_embedded_model)
      .export_values();

  mod.def("parseInput",
          [](const std::string & input_file) {
            _aka::getStaticParser().parse(input_file);
          },
          "Parse an Akantu input file");

  py::class_<_aka::Mesh>(mod, "Mesh")
      .def(py::init<_aka::UInt, const _aka::ID &, const _aka::MemoryID &>(),
           py::arg("spatial_dimension"), py::arg("id") = "mesh",
           py::arg("memory_id") = 0)
      .def("read", &_aka::Mesh::read, py::arg("filename"),
           py::arg("mesh_io_type") = _aka::_miot_auto,
           "read the mesh from a file");

  return mod;
} // namespace

} // namespace
