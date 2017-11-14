/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"
#include "test_gtest_utils.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <vector>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_TEST_SOLID_MECHANICS_MODEL_FIXTURE_HH__
#define __AKANTU_TEST_SOLID_MECHANICS_MODEL_FIXTURE_HH__

using namespace akantu;

template <typename type_> class TestSMMFixture : public ::testing::Test {
public:
  static constexpr const ElementType type = type_::value;
  static constexpr const size_t spatial_dimension = ElementClass<type>::getSpatialDimension();

  void SetUp() override {
    const auto spatial_dimension = this->spatial_dimension;

    std::stringstream element_type;
    element_type << type;
    SCOPED_TRACE(element_type.str().c_str());

    getStaticParser().parse("test_solid_mechanics_model_linear_elastic_"
                            "potential_energy_material.dat");

    mesh = std::make_unique<Mesh>(spatial_dimension);

    mesh->read(element_type.str() + ".msh");

    model = std::make_unique<SolidMechanicsModel>(*mesh, _all_dimensions,
                                                  element_type.str());
  }

  void TearDown() override {
    mesh.reset(nullptr);
    model.reset(nullptr);
  }

protected:
  std::unique_ptr<Mesh> mesh;
  std::unique_ptr<SolidMechanicsModel> model;
};

using types = gtest_list_t<TestElementTypes>;

TYPED_TEST_CASE(TestSMMFixture, types);

#endif /* __AKANTU_TEST_SOLID_MECHANICS_MODEL_FIXTURE_HH__ */
