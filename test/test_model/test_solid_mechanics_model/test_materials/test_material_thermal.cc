/* -------------------------------------------------------------------------- */
#include "material_thermal.hh"
#include "solid_mechanics_model.hh"
#include "test_material_fixtures.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <type_traits>
/* -------------------------------------------------------------------------- */

using namespace akantu;

using types =
    ::testing::Types<Traits<MaterialThermal, 1>, Traits<MaterialThermal, 2>,
                     Traits<MaterialThermal, 3>>;

/* -------------------------------------------------------------------------- */
template <> void FriendMaterial<MaterialThermal<3>>::testComputeStress() {
  Real E = 1.;
  Real nu = .3;
  Real alpha = 2;
  setParam("E", E);
  setParam("nu", nu);
  setParam("alpha", alpha);

  Real deltaT = 1;
  Real sigma = 0;
  this->computeStressOnQuad(sigma, deltaT);
  Real solution = -E / (1 - 2 * nu) * alpha * deltaT;
  auto error = std::abs(sigma - solution);
  ASSERT_NEAR(error, 0, 1e-14);
}

template <> void FriendMaterial<MaterialThermal<2>>::testComputeStress() {
  Real E = 1.;
  Real nu = .3;
  Real alpha = 2;
  setParam("E", E);
  setParam("nu", nu);
  setParam("alpha", alpha);

  Real deltaT = 1;
  Real sigma = 0;
  this->computeStressOnQuad(sigma, deltaT);
  Real solution = -E / (1 - 2 * nu) * alpha * deltaT;
  auto error = std::abs(sigma - solution);
  ASSERT_NEAR(error, 0, 1e-14);
}

template <> void FriendMaterial<MaterialThermal<1>>::testComputeStress() {
  Real E = 1.;
  Real nu = .3;
  Real alpha = 2;
  setParam("E", E);
  setParam("nu", nu);
  setParam("alpha", alpha);

  Real deltaT = 1;
  Real sigma = 0;
  this->computeStressOnQuad(sigma, deltaT);
  Real solution = -E * alpha * deltaT;
  auto error = std::abs(sigma - solution);
  ASSERT_NEAR(error, 0, 1e-14);
}

namespace {

template <typename T>
class TestMaterialThermalFixture : public ::TestMaterialFixture<T> {};

TYPED_TEST_CASE(TestMaterialThermalFixture, types);

TYPED_TEST(TestMaterialThermalFixture, ThermalComputeStress) {
  this->material->testComputeStress();
}
}
