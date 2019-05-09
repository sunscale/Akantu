/* -------------------------------------------------------------------------- */
#include "py_aka_common.hh"
#include "py_aka_boundary_conditions.hh"
#include "py_aka_error.hh"
#include "py_aka_fe_engine.hh"
#include "py_aka_heat_transfer_model.hh"
#include "py_aka_material.hh"
#include "py_aka_mesh.hh"
#include "py_aka_model.hh"
#include "py_aka_parser.hh"
#include "py_aka_solid_mechanics_model.hh"
#include "py_aka_solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */
#include <aka_common.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;

namespace akantu {

/* -------------------------------------------------------------------------- */

void register_initialize(py::module & mod) {
  mod.def("__initialize", []() {
    int nb_args = 0;
    char ** null = nullptr;
    initialize(nb_args, null);
  });
}

/* -------------------------------------------------------------------------- */

void register_enums(py::module & mod) {
  py::enum_<SpatialDirection>(mod, "SpatialDirection")
      .value("_x", _x)
      .value("_y", _y)
      .value("_z", _z)
      .export_values();

  py::enum_<AnalysisMethod>(mod, "AnalysisMethod")
      .value("_static", _static)
      .value("_implicit_dynamic", _implicit_dynamic)
      .value("_explicit_lumped_mass", _explicit_lumped_mass)
      .value("_explicit_lumped_capacity", _explicit_lumped_capacity)
      .value("_explicit_consistent_mass", _explicit_consistent_mass)
      .export_values();

  py::enum_<NonLinearSolverType>(mod, "NonLinearSolverType")
      .value("_nls_linear", _nls_linear)
      .value("_nls_newton_raphson", _nls_newton_raphson)
      .value("_nls_newton_raphson_modified", _nls_newton_raphson_modified)
      .value("_nls_lumped", _nls_lumped)
      .value("_nls_auto", _nls_auto)
      .export_values();

  py::enum_<TimeStepSolverType>(mod, "TimeStepSolverType")
      .value("_tsst_static", _tsst_static)
      .value("_tsst_dynamic", _tsst_dynamic)
      .value("_tsst_dynamic_lumped", _tsst_dynamic_lumped)
      .value("_tsst_not_defined", _tsst_not_defined)
      .export_values();

  py::enum_<IntegrationSchemeType>(mod, "IntegrationSchemeType")
      .value("_ist_pseudo_time", _ist_pseudo_time)
      .value("_ist_forward_euler", _ist_forward_euler)
      .value("_ist_trapezoidal_rule_1", _ist_trapezoidal_rule_1)
      .value("_ist_backward_euler", _ist_backward_euler)
      .value("_ist_central_difference", _ist_central_difference)
      .value("_ist_fox_goodwin", _ist_fox_goodwin)
      .value("_ist_trapezoidal_rule_2", _ist_trapezoidal_rule_2)
      .value("_ist_linear_acceleration", _ist_linear_acceleration)
      .value("_ist_newmark_beta", _ist_newmark_beta)
      .value("_ist_generalized_trapezoidal", _ist_generalized_trapezoidal)
      .export_values();

  py::enum_<SolveConvergenceCriteria>(mod, "SolveConvergenceCriteria")
      .value("_scc_residual", _scc_residual)
      .value("_scc_solution", _scc_solution)
      .value("_scc_residual_mass_wgh", _scc_residual_mass_wgh)
      .export_values();

  py::enum_<CohesiveMethod>(mod, "CohesiveMethod")
      .value("_intrinsic", _intrinsic)
      .value("_extrinsic", _extrinsic)
      .export_values();

  py::enum_<GhostType>(mod, "GhostType")
      .value("_not_ghost", _not_ghost)
      .value("_ghost", _ghost)
      .value("_casper", _casper)
      .export_values();

  py::enum_<MeshIOType>(mod, "MeshIOType")
      .value("_miot_auto", _miot_auto)
      .value("_miot_gmsh", _miot_gmsh)
      .value("_miot_gmsh_struct", _miot_gmsh_struct)
      .value("_miot_diana", _miot_diana)
      .value("_miot_abaqus", _miot_abaqus)
      .export_values();

  py::enum_<ModelType>(mod, "ModelType")
      .value("_model", ModelType::_model)
      .value("_solid_mechanics_model", ModelType::_solid_mechanics_model)
      .value("_solid_mechanics_model_cohesive",
             ModelType::_solid_mechanics_model_cohesive)
      .value("_heat_transfer_model", ModelType::_heat_transfer_model)
      .value("_structural_mechanics_model",
             ModelType::_structural_mechanics_model)
      .value("_embedded_model", ModelType::_embedded_model)
      .export_values();

  py::enum_<ElementType>(mod, "ElementType")
      .value("_point_1", _point_1)
      .value("_segment_2", _segment_2)
      .value("_segment_3", _segment_3)
      .value("_triangle_3", _triangle_3)
      .value("_triangle_6", _triangle_6)
      .value("_quadrangle_4", _quadrangle_4)
      .value("_quadrangle_8", _quadrangle_8)
      .value("_tetrahedron_4", _tetrahedron_4)
      .value("_tetrahedron_10", _tetrahedron_10)
      .value("_pentahedron_6", _pentahedron_6)
      .value("_pentahedron_15", _pentahedron_15)
      .value("_hexahedron_8", _hexahedron_8)
      .value("_hexahedron_20", _hexahedron_20)
      .export_values();

  py::enum_<ElementKind>(mod, "ElementKind")
      .value("_ek_regular", _ek_regular)
      .export_values();

  py::enum_<MatrixType>(mod, "MatrixType")
      .value("_unsymmetric", _unsymmetric)
      .value("_symmetric", _symmetric)
      .export_values();
}

/* -------------------------------------------------------------------------- */
#define AKANTU_PP_STR_TO_TYPE2(s, data, elem) ({BOOST_PP_STRINGIZE(elem), elem})

void register_functions(py::module & mod) {

  mod.def("getElementTypes", []() {
    std::map<std::string, akantu::ElementType> element_types{
        BOOST_PP_SEQ_FOR_EACH_I(
            AKANTU_PP_ENUM, BOOST_PP_SEQ_SIZE(AKANTU_ek_regular_ELEMENT_TYPE),
            BOOST_PP_SEQ_TRANSFORM(AKANTU_PP_STR_TO_TYPE2, akantu,
                                   AKANTU_ek_regular_ELEMENT_TYPE))};

    return element_types;
  });
}

#undef AKANTU_PP_STR_TO_TYPE2

/* -------------------------------------------------------------------------- */

__attribute__((visibility("default"))) void register_all(py::module & mod) {
  register_initialize(mod);
  register_enums(mod);
  register_error(mod);
  register_functions(mod);
  register_parser(mod);
  register_boundary_conditions(mod);
  register_fe_engine(mod);
  register_model(mod);
  register_heat_transfer_model(mod);
  register_solid_mechanics_model(mod);
  register_solid_mechanics_model_cohesive(mod);
  register_material(mod);
  register_mesh(mod);
}

} // namespace akantu
