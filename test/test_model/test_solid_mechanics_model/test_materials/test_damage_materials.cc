/* -------------------------------------------------------------------------- */
#include "material_marigo.hh"
#include "material_mazars.hh"
#include "solid_mechanics_model.hh"
#include "test_material_fixtures.hh"
#include <gtest/gtest.h>
#include <type_traits>
/* -------------------------------------------------------------------------- */

using namespace akantu;

using types = ::testing::Types<

    Traits<MaterialMarigo, 1>, Traits<MaterialMarigo, 2>,
    Traits<MaterialMarigo, 3>,

    Traits<MaterialMazars, 1>, Traits<MaterialMazars, 2>,
    Traits<MaterialMazars, 3>>;


/*****************************************************************/

namespace {

template <typename T>
class TestDamageMaterialFixture : public ::TestMaterialFixture<T> {};

TYPED_TEST_CASE(TestDamageMaterialFixture, types);

TYPED_TEST(TestDamageMaterialFixture, DISABLED_ComputeStress) {
  this->material->testComputeStress();
}
TYPED_TEST(TestDamageMaterialFixture, DISABLED_EnergyDensity) {
  this->material->testEnergyDensity();
}
TYPED_TEST(TestDamageMaterialFixture, DISABLED_ComputeTangentModuli) {
  this->material->testComputeTangentModuli();
}

TYPED_TEST(TestDamageMaterialFixture, DISABLED_ComputePushWaveSpeed) {
  this->material->testPushWaveSpeed();
}

TYPED_TEST(TestDamageMaterialFixture, DISABLED_ComputeShearWaveSpeed) {
  this->material->testShearWaveSpeed();
}
}
/*****************************************************************/
