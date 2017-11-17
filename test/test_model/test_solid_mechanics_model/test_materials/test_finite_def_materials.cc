/* -------------------------------------------------------------------------- */
#include "material_neohookean.hh"
#include "solid_mechanics_model.hh"
#include "test_material_fixtures.hh"
#include <gtest/gtest.h>
#include <type_traits>
/* -------------------------------------------------------------------------- */

using namespace akantu;

using types = ::testing::Types<

    Traits<MaterialNeohookean, 1>, Traits<MaterialNeohookean, 2>,
    Traits<MaterialNeohookean, 3>>;

/*****************************************************************/

template <> void FriendMaterial<MaterialNeohookean<3>>::testComputeStress() {
  TO_IMPLEMENT;
}

/*****************************************************************/
template <>
void FriendMaterial<MaterialNeohookean<3>>::testComputeTangentModuli() {
  TO_IMPLEMENT;
}

/*****************************************************************/

template <> void FriendMaterial<MaterialNeohookean<3>>::testEnergyDensity() {
  TO_IMPLEMENT;
}

/*****************************************************************/

namespace {

template <typename T>
class TestFiniteDefMaterialFixture : public ::TestMaterialFixture<T> {};

TYPED_TEST_CASE(TestFiniteDefMaterialFixture, types);

TYPED_TEST(TestFiniteDefMaterialFixture, FiniteDefComputeStress) {
  this->material->testComputeStress();
}
TYPED_TEST(TestFiniteDefMaterialFixture, FiniteDefEnergyDensity) {
  this->material->testEnergyDensity();
}
TYPED_TEST(TestFiniteDefMaterialFixture, FiniteDefComputeTangentModuli) {
  this->material->testComputeTangentModuli();
}
}
/*****************************************************************/
