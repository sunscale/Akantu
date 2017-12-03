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

template <typename type_> class TestFEMFixture : public ::testing::Test {
public:
  static constexpr const ElementType type = type_::value;
  static constexpr const size_t dim = ElementClass<type>::getSpatialDimension();
  using FEM = FEEngineTemplate<IntegratorGauss, ShapeLagrange>;

  void SetUp() override {
    const auto dim = this->dim;
    const auto type = this->type;
    mesh = std::make_unique<Mesh>(dim);

    std::stringstream meshfilename;
    meshfilename << type << ".msh";
    mesh->read(meshfilename.str());

    lower = mesh->getLowerBounds();
    upper = mesh->getUpperBounds();

    nb_element = this->mesh->getNbElement(type);

    fem = std::make_unique<FEM>(*mesh, dim, "my_fem");
    fem->initShapeFunctions();

    nb_quadrature_points_total = fem->getNbIntegrationPoints(type) * nb_element;

    std::stringstream sstr;
    sstr << type;
    SCOPED_TRACE(sstr.str().c_str());
  }

  void TearDown() override {
    fem.reset(nullptr);
    mesh.reset(nullptr);
  }

protected:
  std::unique_ptr<FEM> fem;
  std::unique_ptr<Mesh> mesh;
  UInt nb_element;
  UInt nb_quadrature_points_total;
  Vector<Real> lower;
  Vector<Real> upper;
};

template <typename type_>
constexpr const ElementType TestFEMFixture<type_>::type;

template <typename type_>
constexpr const size_t TestFEMFixture<type_>::dim;


using types = gtest_list_t<TestElementTypes>;

TYPED_TEST_CASE(TestFEMFixture, types);

#endif /* __AKANTU_TEST_FE_ENGINE_FIXTURE_HH__ */
