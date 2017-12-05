/* -------------------------------------------------------------------------- */
#include "mesh_utils.hh"
#include "non_linear_solver.hh"
#include "solid_mechanics_model.hh"
#include "test_gtest_utils.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <vector>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_PATCH_TEST_LINEAR_FIXTURE_HH__
#define __AKANTU_PATCH_TEST_LINEAR_FIXTURE_HH__

//#define DEBUG_TEST

using namespace akantu;

template <typename tuple_>
class TestPatchTestLinear : public ::testing::Test {
public:
  static constexpr ElementType type = std::tuple_element_t<0, tuple_>::value;
  static constexpr bool plane_strain = std::tuple_element_t<1, tuple_>::value;
  static constexpr size_t dim = ElementClass<type>::getSpatialDimension();

  void SetUp() override {
    std::stringstream element_type;
    element_type << type;

    mesh = std::make_unique<Mesh>(dim);
    mesh->read(element_type.str() + ".msh");
    MeshUtils::buildFacets(*mesh);
    mesh->createBoundaryGroupFromGeometry();

    if (dim == 2)
      element_type << ":plane_" << (plane_strain ? "stain" : "stress");

    SCOPED_TRACE(element_type.str().c_str());

    model = std::make_unique<SolidMechanicsModel>(*mesh);

#ifdef DEBUG_TEST
    model->solveStep();
    model->addDumpField("strain");
    model->addDumpField("stress");
    model->addDumpField("external_force");
    model->addDumpField("internal_force");
    model->addDumpField("velocity");
    model->addDumpField("acceleration");
    model->addDumpField("displacement");
    model->dump();
#endif
  }

  void TearDown() override {
    model.reset(nullptr);
    mesh.reset(nullptr);
  }

  void initModel(const AnalysisMethod & method, const std::string & material_file) {
    debug::setDebugLevel(dblError);
    getStaticParser().parse(material_file);

    this->model->initFull(_analysis_method = method);

    this->applyBC();
  }

  void applyBC() {
    const auto & coordinates = mesh->getNodes();
    auto & displacement = model->getDisplacement();
    auto & boundary = model->getBlockedDOFs();

    for (auto & eg : mesh->getElementGroups())
      for (const auto & node : eg.second->getNodeGroup()) {
        for (UInt s = 0; s < dim; ++s) {
          boundary(node, s) = true;
        }

        setDisplacement(displacement.begin(dim)[node], coordinates.begin(dim)[node]);
      }
  }

  template<typename V1, typename V2>
  void setDisplacement(V1 && disp, V2 && coord) {
    for (UInt i = 0; i < dim; ++i) {
      disp(i) = this->alpha(i, 0);
      for (UInt j = 0; j < dim; ++j) {
        disp(i) += this->alpha(i, j + 1) * coord(j);
      }
    }
  }

  Matrix<Real> prescribed_strain() {
    Matrix<Real> strain(dim, dim);

    for (UInt i = 0; i < dim; ++i) {
      for (UInt j = 0; j < dim; ++j) {
        strain(i, j) = alpha(i, j + 1);
      }
    }
    return strain;
  }

  void checkStrains() {
    auto pstrain = prescribed_strain();

    for (auto & strain :
             make_view(this->model->getMaterial(0).getGradU(type), dim, dim)) {
      auto diff = strain - pstrain;
      auto strain_error = diff.template norm<L_inf>() / strain.template norm<L_inf>();

      EXPECT_NEAR(0, strain_error, strain_tolerance);
    }
  }

  template<typename pstress_func_t>
  void checkStresses(pstress_func_t && pstress_func) {
    auto pstress = pstress_func(prescribed_strain());

    for (auto & stress :
             make_view(this->model->getMaterial(0).getStress(type), dim, dim)) {
      auto diff = stress - pstress;
      auto stress_error = diff.template norm<L_inf>() / stress.template norm<L_inf>();

      EXPECT_NEAR(0, stress_error, stress_tolerance);
    }
  }

  void checkDisplacements() {
    const auto & coordinates = mesh->getNodes();
    auto & displacement = model->getDisplacement();
    Vector<Real> ref_disp(dim);

    for (auto && tuple :
         zip(make_view(coordinates, dim), make_view(displacement, dim))) {
      setDisplacement(ref_disp, std::get<0>(tuple));
      auto diff = std::get<1>(tuple) - ref_disp;
      auto disp_error = diff.template norm<L_inf>();

      EXPECT_NEAR(0, disp_error, disp_tolerance);
    }
  }

protected:
  std::unique_ptr<Mesh> mesh;
  std::unique_ptr<SolidMechanicsModel> model;
  Matrix<Real> alpha{{0.01, 0.02, 0.03, 0.04},
                     {0.05, 0.06, 0.07, 0.08},
                     {0.09, 0.10, 0.11, 0.12}};

  Real strain_tolerance{1e-14};
  Real stress_tolerance{1e-14};
  Real disp_tolerance{1e-15};
};

template <typename tuple_>
constexpr ElementType TestPatchTestLinear<tuple_>::type;

template <typename tuple_>
constexpr bool TestPatchTestLinear<tuple_>::plane_strain;

template <typename tuple_>
constexpr size_t TestPatchTestLinear<tuple_>::dim;

template <typename T> struct invalid_plan_stress : std::true_type {};

template <typename type, typename bool_c>
struct invalid_plan_stress<std::tuple<type, bool_c>>
    : aka::bool_constant<ElementClass<type::value>::getSpatialDimension() !=
                             2 and
                         not bool_c::value> {};

using true_false =
    std::tuple<aka::bool_constant<true>, aka::bool_constant<false>>;

template <typename T> using valid_types = aka::negation<invalid_plan_stress<T>>;

using types = gtest_list_t<
    tuple_filter_t<valid_types, cross_product_t<TestElementTypes, true_false>>>;

TYPED_TEST_CASE(TestPatchTestLinear, types);

#endif /* __AKANTU_PATCH_TEST_LINEAR_FIXTURE_HH__ */
