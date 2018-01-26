/* -------------------------------------------------------------------------- */
#include "element_class.hh"
#include "fe_engine.hh"
#include "integrator_gauss.hh"
#include "shape_lagrange.hh"
#include "test_gtest_utils.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_TEST_FE_ENGINE_FIXTURE_HH__
#define __AKANTU_TEST_FE_ENGINE_FIXTURE_HH__

using namespace akantu;

/// Generic class for FEEngine tests
template <typename type_, template <ElementKind> class shape_t,
          ElementKind kind = _ek_regular>
class TestFEMBaseFixture : public ::testing::Test {
public:
  static constexpr const ElementType type = type_::value;
  static constexpr const size_t dim = ElementClass<type>::getSpatialDimension();
  using FEM = FEEngineTemplate<IntegratorGauss, shape_t, kind>;

  /// Setup reads mesh corresponding to element type and initializes an FEEngine
  void SetUp() override {
    const auto dim = this->dim;
    const auto type = this->type;
    mesh = std::make_unique<Mesh>(dim);

    std::stringstream meshfilename;
    meshfilename << type << ".msh";
    this->readMesh(meshfilename.str());

    lower = mesh->getLowerBounds();
    upper = mesh->getUpperBounds();

    nb_element = this->mesh->getNbElement(type);

    fem = std::make_unique<FEM>(*mesh, dim, "my_fem");
    fem->initShapeFunctions();

    nb_quadrature_points_total = fem->getNbIntegrationPoints(type) * nb_element;

    SCOPED_TRACE(aka::to_string(type));
  }

  void TearDown() override {
    fem.reset(nullptr);
    mesh.reset(nullptr);
  }

  /// Should be reimplemented if further treatment of the mesh is needed
  virtual void readMesh(std::string file_name) { mesh->read(file_name); }

protected:
  std::unique_ptr<FEM> fem;
  std::unique_ptr<Mesh> mesh;
  UInt nb_element;
  UInt nb_quadrature_points_total;
  Vector<Real> lower;
  Vector<Real> upper;
};

template <typename type_, template <ElementKind> class shape_t,
          ElementKind kind>
constexpr const ElementType TestFEMBaseFixture<type_, shape_t, kind>::type;

template <typename type_, template <ElementKind> class shape_t,
          ElementKind kind>
constexpr const size_t TestFEMBaseFixture<type_, shape_t, kind>::dim;

/* -------------------------------------------------------------------------- */
/// Base class for test with Lagrange FEEngine and regular elements
template <typename type_>
using TestFEMFixture = TestFEMBaseFixture<type_, ShapeLagrange>;

/* -------------------------------------------------------------------------- */

using types = gtest_list_t<TestElementTypes>;

TYPED_TEST_CASE(TestFEMFixture, types);

#endif /* __AKANTU_TEST_FE_ENGINE_FIXTURE_HH__ */
