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
    Traits<MaterialElasticOrthotropic, 3>,

    Traits<MaterialElasticLinearAnisotropic, 1>,
    Traits<MaterialElasticLinearAnisotropic, 2>,
    Traits<MaterialElasticLinearAnisotropic, 3>>;

/*****************************************************************/

template <> void FriendMaterial<MaterialElastic<3>>::testComputeStress() {

  Real E = 1.;
  Real nu = .3;
  setParam("E", E);
  setParam("nu", nu);
  setParam("rho", nu);

  Matrix<Real> rotation_matrix = getRandomRotation3d();

  auto grad_u = this->getDeviatoricStrain(1.);

  auto grad_u_rot = this->applyRotation(grad_u, rotation_matrix);

  Matrix<Real> sigma_rot(3, 3);
  this->computeStressOnQuad(grad_u_rot, sigma_rot, 0.);

  auto sigma = this->reverseRotation(sigma_rot, rotation_matrix);

  Matrix<Real> sigma_expected =
    0.5 * E / (1 + nu) * (grad_u + grad_u.transpose());

  auto diff = sigma - sigma_expected;
  Real stress_error = diff.norm<L_inf>();
  ASSERT_NEAR(stress_error, 0., 1e-14);
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
