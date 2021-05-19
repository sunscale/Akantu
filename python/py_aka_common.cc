/* -------------------------------------------------------------------------- */
#include <aka_common.hh>
/* -------------------------------------------------------------------------- */
#include <boost/preprocessor.hpp>
/* -------------------------------------------------------------------------- */
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;

namespace akantu {

/* -------------------------------------------------------------------------- */
#define PY_AKANTU_PP_VALUE(s, data, elem)                                      \
  .value(BOOST_PP_STRINGIZE(elem), BOOST_PP_CAT(data, elem))

#define PY_AKANTU_REGISTER_ENUM_(type_name, list, prefix, mod)                 \
  py::enum_<type_name>(mod, BOOST_PP_STRINGIZE(type_name))                     \
      BOOST_PP_SEQ_FOR_EACH(PY_AKANTU_PP_VALUE, prefix, list)                  \
          .export_values()

#define PY_AKANTU_REGISTER_CLASS_ENUM(type_name, list, mod)                    \
  PY_AKANTU_REGISTER_ENUM_(type_name, list, type_name::_, mod)

#define PY_AKANTU_REGISTER_ENUM(type_name, list, mod)                          \
  PY_AKANTU_REGISTER_ENUM_(type_name, list, , mod)

/* -------------------------------------------------------------------------- */
void register_initialize(py::module & mod) {
  mod.def("__initialize", []() {
    int nb_args = 0;
    char ** null = nullptr;
    initialize(nb_args, null);
  });
}

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
  

  PY_AKANTU_REGISTER_CLASS_ENUM(ModelType, AKANTU_MODEL_TYPES, mod);
  PY_AKANTU_REGISTER_CLASS_ENUM(NonLinearSolverType,
                                AKANTU_NON_LINEAR_SOLVER_TYPES, mod);
  PY_AKANTU_REGISTER_CLASS_ENUM(TimeStepSolverType,
                                AKANTU_TIME_STEP_SOLVER_TYPE, mod);
  PY_AKANTU_REGISTER_CLASS_ENUM(IntegrationSchemeType,
                                AKANTU_INTEGRATION_SCHEME_TYPE, mod);
  PY_AKANTU_REGISTER_CLASS_ENUM(SolveConvergenceCriteria,
                                AKANTU_SOLVE_CONVERGENCE_CRITERIA, mod);

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

  py::enum_<MatrixType>(mod, "MatrixType")
      .value("_unsymmetric", _unsymmetric)
      .value("_symmetric", _symmetric)
      .export_values();

  PY_AKANTU_REGISTER_ENUM(ElementType, AKANTU_ALL_ELEMENT_TYPE(_not_defined),
                          mod);
  PY_AKANTU_REGISTER_ENUM(ElementKind, AKANTU_ELEMENT_KIND(_ek_not_defined),
                          mod);
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

} // namespace akantu
