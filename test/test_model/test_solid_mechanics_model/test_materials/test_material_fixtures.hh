#include "material_elastic.hh"
#include "solid_mechanics_model.hh"
#include <gtest/gtest.h>
#include <type_traits>
/* -------------------------------------------------------------------------- */

#define TO_IMPLEMENT AKANTU_EXCEPTION("TEST TO IMPLEMENT")

using namespace akantu;

/*****************************************************************/

template <typename T> class FriendMaterial : public T {

public:
  virtual void testComputeStress() { TO_IMPLEMENT; };
  virtual void testComputeTangentModuli() { TO_IMPLEMENT; };
  virtual void testEnergyDensity() { TO_IMPLEMENT; };
  virtual void testPushWaveSpeed() { TO_IMPLEMENT; }
  virtual void testShearWaveSpeed() { TO_IMPLEMENT; }  

  FriendMaterial(SolidMechanicsModel & model, const ID & id = "")
      : T(model, id) {}
};

/*****************************************************************/

template <typename T> class TestMaterialFixture : public ::testing::Test {

public:
  using mat_class = typename T::mat_class;

  void SetUp() override {
    constexpr auto spatial_dimension = T::Dim;
    mesh = std::make_unique<Mesh>(spatial_dimension);
    model = std::make_unique<SolidMechanicsModel>(*mesh);
    material = std::make_unique<friend_class>(*model);
  }

  void TearDown() override {
    material.reset(nullptr);
    model.reset(nullptr);
    mesh.reset(nullptr);
  }

  using friend_class = FriendMaterial<mat_class>;

protected:
  std::unique_ptr<Mesh> mesh;
  std::unique_ptr<SolidMechanicsModel> model;
  std::unique_ptr<friend_class> material;
};

/*****************************************************************/
template <template <UInt> class T, UInt _Dim> struct Traits {
  using mat_class = T<_Dim>;
  static constexpr UInt Dim = _Dim;
};
/*****************************************************************/
