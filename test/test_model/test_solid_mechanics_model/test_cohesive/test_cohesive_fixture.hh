/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model_cohesive.hh"
#include "test_gtest_utils.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <vector>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_TEST_COHESIVE_FIXTURE_HH__
#define __AKANTU_TEST_COHESIVE_FIXTURE_HH__

using namespace akantu;

template <typename param_> class TestSMMCFixture : public ::testing::Test {
public:
  static constexpr ElementType cohesive_type =
      std::tuple_element_t<0, param_>::value;
  static constexpr size_t dim =
      ElementClass<cohesive_type>::getSpatialDimension();
  static constexpr bool is_extrinsic = std::tuple_element_t<1, param_>::value;

  void SetUp() override {
    getStaticParser().parse(this->input_file);

    mesh = std::make_unique<Mesh>(this->dim);
    EXPECT_NO_THROW({ mesh->read(this->mesh_name); });

    model = std::make_unique<SolidMechanicsModelCohesive>(*mesh);
    model->initFull(_analysis_method = this->analysis_method,
                    _is_extrinsic = this->is_extrinsic);

    model->applyBC(BC::Dirichlet::FlagOnly(_x), "fixed");
    if (this->dim > 1)
      model->applyBC(BC::Dirichlet::FlagOnly(_y), "fixed");
    if (this->dim > 2)
      model->applyBC(BC::Dirichlet::FlagOnly(_z), "fixed");

    auto & fe_engine = model->getFEEngine("FacetsFEEngine");
    auto & mesh_facets = fe_engine.getMesh();
    auto facet_type = mesh_facets.getFacetType(this->cohesive_type);
    auto & group =
        mesh_facets.getElementGroup("insertion").getElements(facet_type);
    Array<Real> ones(fe_engine.getNbIntegrationPoints(facet_type) *
                     group.size());
    ones.set(1.);

    surface = fe_engine.integrate(ones, facet_type, _not_ghost, group);
  }

  void TearDown() override {
    model.reset(nullptr);
    mesh.reset(nullptr);
  }

protected:
  std::unique_ptr<Mesh> mesh;
  std::unique_ptr<SolidMechanicsModelCohesive> model;
  AnalysisMethod analysis_method{_explicit_lumped_mass};
  std::string mesh_name{aka::to_string(cohesive_type) + ".msh"};
  std::string input_file{"material.dat"};

  Real surface{0};
};

/* -------------------------------------------------------------------------- */
template <typename param_>
constexpr ElementType TestSMMCFixture<param_>::cohesive_type;
template <typename param_> constexpr size_t TestSMMCFixture<param_>::dim;
template <typename param_> constexpr bool TestSMMCFixture<param_>::is_extrinsic;
/* -------------------------------------------------------------------------- */

using IsExtrinsicTypes = std::tuple<std::true_type, std::false_type>;
using types =
    gtest_list_t<cross_product_t<TestCohesiveElementTypes, IsExtrinsicTypes>>;

TYPED_TEST_CASE(TestSMMCFixture, types);

#endif /* __AKANTU_TEST_COHESIVE_FIXTURE_HH__ */
