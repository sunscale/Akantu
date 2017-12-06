/* -------------------------------------------------------------------------- */
#include "structural_mechanics_model.hh"
#include "element_class_structural.hh"
#include "test_gtest_utils.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <vector>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_TEST_STRUCTURAL_MECHANICS_MODEL_FIXTURE_HH__
#define __AKANTU_TEST_STRUCTURAL_MECHANICS_MODEL_FIXTURE_HH__

using namespace akantu;

template <typename type_> class TestStructuralFixture : public ::testing::Test {
public:
  static constexpr const ElementType type = type_::value;
  static constexpr const size_t spatial_dimension =
      ElementClass<type>::getSpatialDimension();

  void SetUp() override {
    const auto spatial_dimension = this->spatial_dimension;

    mesh = std::make_unique<Mesh>(spatial_dimension);

    // mesh->read(this->makeMeshName());

    std::stringstream element_type;
    element_type << type;

    model = std::make_unique<StructuralMechanicsModel>(*mesh, _all_dimensions,
                                                       element_type.str());
  }

  virtual std::string makeMeshName() {
    std::stringstream element_type;
    element_type << type;
    SCOPED_TRACE(element_type.str().c_str());
    return element_type.str() + ".msh";
  }

  void TearDown() override {
    model.reset(nullptr);
    mesh.reset(nullptr);
  }

protected:
  std::unique_ptr<Mesh> mesh;
  std::unique_ptr<StructuralMechanicsModel> model;
};

using types = gtest_list_t<StructuralTestElementTypes>;

TYPED_TEST_CASE(TestStructuralFixture, types);

#endif /* __AKANTU_TEST_STRUCTURAL_MECHANICS_MODEL_FIXTURE_HH__ */
