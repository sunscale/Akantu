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

// This fixture uses very small meshes with a volume of 1.
template <typename type_> class TestSMMFixture : public ::testing::Test {
public:
  static constexpr ElementType type = type_::value;
  static constexpr size_t spatial_dimension =
      ElementClass<type>::getSpatialDimension();

  void SetUp() override {
    mesh = std::make_unique<Mesh>(this->spatial_dimension);

    mesh->read(this->mesh_file);

    SCOPED_TRACE(aka::to_string(this->type).c_str());

    model = std::make_unique<SolidMechanicsModel>(*mesh, _all_dimensions,
                                                  aka::to_string(this->type));
  }

  void TearDown() override {
    model.reset(nullptr);
    mesh.reset(nullptr);
  }

protected:
  std::string mesh_file{aka::to_string(this->type) + ".msh"};
  std::unique_ptr<Mesh> mesh;
  std::unique_ptr<SolidMechanicsModel> model;
};

template <typename type_> constexpr ElementType TestSMMFixture<type_>::type;
template <typename type_> constexpr size_t TestSMMFixture<type_>::spatial_dimension;

using types = gtest_list_t<TestElementTypes>;

TYPED_TEST_CASE(TestSMMFixture, types);


// This fixture uses somewhat finer meshes representing bars.
template <typename type_>
class TestSMMFixtureBar : public TestSMMFixture<type_> {
public:
  void SetUp() override {
    this->mesh_file = "bar" + aka::to_string(this->type) + ".msh";
    TestSMMFixture<type_>::SetUp();
  }
};

TYPED_TEST_CASE(TestSMMFixtureBar, types);

#endif /* __AKANTU_TEST_SOLID_MECHANICS_MODEL_FIXTURE_HH__ */
