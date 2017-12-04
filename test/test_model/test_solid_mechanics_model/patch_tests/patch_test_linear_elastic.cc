/* -------------------------------------------------------------------------- */
#include "mesh_utils.hh"
#include "non_linear_solver.hh"
#include "solid_mechanics_model.hh"
#include "test_gtest_utils.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <vector>
/* -------------------------------------------------------------------------- */

//#define DEBUG_TEST

using namespace akantu;

template <typename tuple_>
class TestPatchTestLinearElastic : public ::testing::Test {
public:
  static constexpr ElementType type = std::tuple_element_t<0, tuple_>::value;
  static constexpr bool plane_strain = std::tuple_element_t<1, tuple_>::value;
  static constexpr size_t dim = ElementClass<type>::getSpatialDimension();

  void SetUp() override {
    debug::setDebugLevel(dblError);
    if (plane_strain)
      getStaticParser().parse("material_check_stress_plane_strain.dat");
    else
      getStaticParser().parse("material_check_stress_plane_stress.dat");

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

  Matrix<Real> prescribed_stress() {
    Matrix<Real> stress(dim, dim);

    // plane strain in 2d
    Matrix<Real> strain(dim, dim);
    Matrix<Real> pstrain = prescribed_strain();

    Real nu = this->model->getMaterial(0).get("nu");
    Real E = this->model->getMaterial(0).get("E");

    strain = (pstrain + pstrain.transpose()) / 2.;
    Real trace = strain.trace();

    Real lambda = nu * E / ((1 + nu) * (1 - 2 * nu));
    Real mu = E / (2 * (1 + nu));

    if (not plane_strain) {
      lambda = nu * E / (1 - nu * nu);
    }

    if (dim == 1) {
      stress(0, 0) = E * strain(0, 0);
    } else {
      for (UInt i = 0; i < dim; ++i)
        for (UInt j = 0; j < dim; ++j)
          stress(i, j) = (i == j) * lambda * trace + 2 * mu * strain(i, j);
    }

    return stress;
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

  void checkStresses() {
    auto pstress = prescribed_stress();

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
constexpr ElementType TestPatchTestLinearElastic<tuple_>::type;

template <typename tuple_>
constexpr bool TestPatchTestLinearElastic<tuple_>::plane_strain;

template <typename tuple_>
constexpr size_t TestPatchTestLinearElastic<tuple_>::dim;

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

TYPED_TEST_CASE(TestPatchTestLinearElastic, types);

TYPED_TEST(TestPatchTestLinearElastic, Implicit) {
  this->model->initFull(_analysis_method = _static);

  this->applyBC();

  auto & solver = this->model->getNonLinearSolver();
  solver.set("max_iterations", 2);
  solver.set("threshold", 2e-4);
  solver.set("convergence_type", _scc_residual);

  this->model->solveStep();

  this->checkDisplacements();
  this->checkStrains();
  this->checkStresses();
}

TYPED_TEST(TestPatchTestLinearElastic, Explicit) {
  this->model->initFull(_analysis_method = _explicit_lumped_mass);

  this->applyBC();

  const auto & coordinates = this->mesh->getNodes();
  auto & displacement = this->model->getDisplacement();

  //set the position of all nodes to the static solution
  for (auto && tuple :
           zip(make_view(coordinates, this->dim), make_view(displacement, this->dim))) {
    this->setDisplacement(std::get<1>(tuple), std::get<0>(tuple));
  }

  for(UInt s = 0; s < 100; ++s) {
    this->model->solveStep();
  }

  auto ekin = this->model->getEnergy("kinetic");
  EXPECT_NEAR(0, ekin, 1e-16);

  this->checkDisplacements();
  this->checkStrains();
  this->checkStresses();
}
