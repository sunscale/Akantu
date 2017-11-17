/* -------------------------------------------------------------------------- */
#include "material_elastic.hh"
#include "material_elastic_orthotropic.hh"
#include "solid_mechanics_model.hh"
#include "test_material_fixtures.hh"
#include <gtest/gtest.h>
#include <type_traits>
/* -------------------------------------------------------------------------- */

using namespace akantu;

using types = ::testing::Types<

    Traits<MaterialElastic, 1>, Traits<MaterialElastic, 2>,
    Traits<MaterialElastic, 3>,

    Traits<MaterialElasticOrthotropic, 1>,
    Traits<MaterialElasticOrthotropic, 2>,
    Traits<MaterialElasticOrthotropic, 3>>;

/*****************************************************************/

template <> void FriendMaterial<MaterialElastic<3>>::testComputeStress() {

  Real E = 1.;
  Real nu = .3;
  setParam("E", E);
  setParam("nu", nu);
  setParam("rho", nu);

  Real epsilon = 1.;
  const Matrix<Real> grad_u = {{0, epsilon, 0}, {epsilon, 0, 0}, {0, 0, 0}};
  Matrix<Real> sigma(3, 3);
  this->computeStressOnQuad(grad_u, sigma, 0.);

  Matrix<Real> sigma_expected = {{0, 1, 0}, {1, 0, 0}, {0, 0, 0}};
  sigma_expected *= E / (1 + nu);

  auto diff = sigma - sigma_expected;
  Real stress_error = diff.norm<L_inf>();
  ASSERT_DOUBLE_EQ(stress_error, 0.);
}

/*****************************************************************/
template <>
void FriendMaterial<MaterialElastic<3>>::testComputeTangentModuli() {
  TO_IMPLEMENT;
}

/*****************************************************************/

template <> void FriendMaterial<MaterialElastic<3>>::testEnergyDensity() {
  TO_IMPLEMENT;
}

/*****************************************************************/

template <> void FriendMaterial<MaterialElastic<2>>::testComputeStress() {
  TO_IMPLEMENT;
}

/*****************************************************************/
template <>
void FriendMaterial<MaterialElastic<2>>::testComputeTangentModuli() {
  TO_IMPLEMENT;
}

/*****************************************************************/

template <> void FriendMaterial<MaterialElastic<2>>::testEnergyDensity() {
  TO_IMPLEMENT;
}

/*****************************************************************/

template <> void FriendMaterial<MaterialElastic<1>>::testComputeStress() {
  TO_IMPLEMENT;
}

/*****************************************************************/
template <>
void FriendMaterial<MaterialElastic<1>>::testComputeTangentModuli() {
  TO_IMPLEMENT;
}

/*****************************************************************/

template <> void FriendMaterial<MaterialElastic<1>>::testEnergyDensity() {
  TO_IMPLEMENT;
}

/*****************************************************************/

namespace {

template <typename T>
class TestElasticMaterialFixture : public ::TestMaterialFixture<T> {};

TYPED_TEST_CASE(TestElasticMaterialFixture, types);

TYPED_TEST(TestElasticMaterialFixture, ElasticComputeStress) {
  this->material->testComputeStress();
}
TYPED_TEST(TestElasticMaterialFixture, ElasticEnergyDensity) {
  this->material->testEnergyDensity();
}
TYPED_TEST(TestElasticMaterialFixture, ElasticComputeTangentModuli) {
  this->material->testComputeTangentModuli();
}

TYPED_TEST(TestElasticMaterialFixture, ElasticComputePushWaveSpeed) {
  this->material->testPushWaveSpeed();
}

TYPED_TEST(TestElasticMaterialFixture, ElasticComputeShearWaveSpeed) {
  this->material->testShearWaveSpeed();
}
}

/*****************************************************************/
