/* -------------------------------------------------------------------------- */
#include "material_elastic.hh"
#include "material_elastic_orthotropic.hh"
#include "solid_mechanics_model.hh"
#include "test_material_fixtures.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <type_traits>
/* -------------------------------------------------------------------------- */

using namespace akantu;

using types =
    ::testing::Types<Traits<MaterialElastic, 1>, Traits<MaterialElastic, 2>,
                     Traits<MaterialElastic, 3>,

                     Traits<MaterialElasticOrthotropic, 1>,
                     Traits<MaterialElasticOrthotropic, 2>,
                     Traits<MaterialElasticOrthotropic, 3>,

                     Traits<MaterialElasticLinearAnisotropic, 1>,
                     Traits<MaterialElasticLinearAnisotropic, 2>,
                     Traits<MaterialElasticLinearAnisotropic, 3>>;

/* -------------------------------------------------------------------------- */
template <> void FriendMaterial<MaterialElastic<3>>::testPushWaveSpeed() {
  Real E = 1.;
  Real nu = .3;
  Real rho = 2;
  setParam("E", E);
  setParam("nu", nu);
  setParam("rho", rho);
  auto wave_speed = this->getPushWaveSpeed(Element());

  Real K = E / (3 * (1 - 2 * nu));
  Real mu = E / (2 * (1 + nu));
  Real sol = std::sqrt((K + 4. / 3 * mu) / rho);

  ASSERT_NEAR(wave_speed, sol, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <> void FriendMaterial<MaterialElastic<3>>::testShearWaveSpeed() {
  Real E = 1.;
  Real nu = .3;
  Real rho = 2;
  setParam("E", E);
  setParam("nu", nu);
  setParam("rho", rho);
  auto wave_speed = this->getShearWaveSpeed(Element());

  Real mu = E / (2 * (1 + nu));
  Real sol = std::sqrt(mu / rho);

  ASSERT_NEAR(wave_speed, sol, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <> void FriendMaterial<MaterialElastic<3>>::testComputeStress() {
  Real E = 1.;
  Real nu = .3;
  Real sigma_th = 0.3; // thermal stress
  setParam("E", E);
  setParam("nu", nu);

  Matrix<Real> rotation_matrix = getRandomRotation3d();

  auto grad_u = this->getDeviatoricStrain(1.);

  auto grad_u_rot = this->applyRotation(grad_u, rotation_matrix);

  Matrix<Real> sigma_rot(3, 3);
  this->computeStressOnQuad(grad_u_rot, sigma_rot, sigma_th);

  auto sigma = this->reverseRotation(sigma_rot, rotation_matrix);

  Matrix<Real> identity(3, 3);
  identity.eye();

  Matrix<Real> sigma_expected =
      0.5 * E / (1 + nu) * (grad_u + grad_u.transpose()) + sigma_th * identity;

  auto diff = sigma - sigma_expected;
  Real stress_error = diff.norm<L_inf>();
  ASSERT_NEAR(stress_error, 0., 1e-14);
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<MaterialElastic<3>>::testComputeTangentModuli() {
  Real E = 1.;
  Real nu = .3;
  setParam("E", E);
  setParam("nu", nu);
  Matrix<Real> tangent(6, 6);

  // clang-format off
  Matrix<Real> solution = {
      {1 - nu, nu, nu, 0, 0, 0},
      {nu, 1 - nu, nu, 0, 0, 0},
      {nu, nu, 1 - nu, 0, 0, 0},
      {0, 0, 0, (1 - 2 * nu) / 2, 0, 0},
      {0, 0, 0, 0, (1 - 2 * nu) / 2, 0},
      {0, 0, 0, 0, 0, (1 - 2 * nu) / 2},
  };
  // clang-format on
  solution *= E / ((1 + nu) * (1 - 2 * nu));

  this->computeTangentModuliOnQuad(tangent);
  Real tangent_error = (tangent - solution).norm<L_2>();
  ASSERT_NEAR(tangent_error, 0, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <> void FriendMaterial<MaterialElastic<3>>::testEnergyDensity() {
  Matrix<Real> sigma = {{1, 2, 3}, {2, 4, 5}, {3, 5, 6}};
  Matrix<Real> eps = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  Real epot = 0;
  Real solution = 5.5;
  this->computePotentialEnergyOnQuad(eps, sigma, epot);
  ASSERT_NEAR(epot, solution, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <> void FriendMaterial<MaterialElastic<1>>::testComputeStress() {
  Real E = 3.;
  setParam("E", E);

  Matrix<Real> eps = {{2}};
  Matrix<Real> sigma(1, 1);
  Real sigma_th = 2;
  this->computeStressOnQuad(eps, sigma, sigma_th);

  auto solution = E * eps(0, 0) + sigma_th;
  ASSERT_NEAR(sigma(0, 0), solution, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <> void FriendMaterial<MaterialElastic<1>>::testEnergyDensity() {
  Real eps = 2, sigma = 2;
  Real epot = 0;
  this->computePotentialEnergyOnQuad({{eps}}, {{sigma}}, epot);
  Real solution = 2;
  ASSERT_NEAR(epot, solution, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<MaterialElastic<1>>::testComputeTangentModuli() {
  Real E = 2;
  setParam("E", E);
  Matrix<Real> tangent(1, 1);
  this->computeTangentModuliOnQuad(tangent);
  ASSERT_NEAR(tangent(0, 0), E, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <> void FriendMaterial<MaterialElastic<1>>::testPushWaveSpeed() {
  Real E = 3., rho = 2.;
  setParam("E", E);
  setParam("rho", rho);

  auto wave_speed = this->getPushWaveSpeed(Element());
  auto solution = std::sqrt(E / rho);
  ASSERT_NEAR(wave_speed, solution, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <> void FriendMaterial<MaterialElastic<1>>::testShearWaveSpeed() {
  ASSERT_THROW(this->getShearWaveSpeed(Element()), debug::Exception);
}

/* -------------------------------------------------------------------------- */
template <> void FriendMaterial<MaterialElastic<2>>::testPushWaveSpeed() {
  Real E = 1.;
  Real nu = .3;
  Real rho = 2;
  setParam("E", E);
  setParam("nu", nu);
  setParam("rho", rho);
  auto wave_speed = this->getPushWaveSpeed(Element());

  Real K = E / (3 * (1 - 2 * nu));
  Real mu = E / (2 * (1 + nu));
  Real sol = std::sqrt((K + 4. / 3 * mu) / rho);

  ASSERT_NEAR(wave_speed, sol, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <> void FriendMaterial<MaterialElastic<2>>::testShearWaveSpeed() {
  Real E = 1.;
  Real nu = .3;
  Real rho = 2;
  setParam("E", E);
  setParam("nu", nu);
  setParam("rho", rho);
  auto wave_speed = this->getShearWaveSpeed(Element());

  Real mu = E / (2 * (1 + nu));
  Real sol = std::sqrt(mu / rho);

  ASSERT_NEAR(wave_speed, sol, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <> void FriendMaterial<MaterialElastic<2>>::testComputeStress() {
  Real E = 1.;
  Real nu = .3;
  Real sigma_th = 0.3; // thermal stress
  setParam("E", E);
  setParam("nu", nu);

  Matrix<Real> rotation_matrix = getRandomRotation2d();

  auto grad_u = this->getDeviatoricStrain(1.).block(0, 0, 2, 2);

  auto grad_u_rot = this->applyRotation(grad_u, rotation_matrix);

  Matrix<Real> sigma_rot(2, 2);
  this->computeStressOnQuad(grad_u_rot, sigma_rot, sigma_th);

  auto sigma = this->reverseRotation(sigma_rot, rotation_matrix);

  Matrix<Real> identity(2, 2);
  identity.eye();

  Matrix<Real> sigma_expected =
      0.5 * E / (1 + nu) * (grad_u + grad_u.transpose()) + sigma_th * identity;

  auto diff = sigma - sigma_expected;
  Real stress_error = diff.norm<L_inf>();
  ASSERT_NEAR(stress_error, 0., 1e-14);
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<MaterialElastic<2>>::testComputeTangentModuli() {
  Real E = 1.;
  Real nu = .3;
  setParam("E", E);
  setParam("nu", nu);
  Matrix<Real> tangent(3, 3);

  /* Plane Strain */
  // clang-format off
  Matrix<Real> solution = {
    {1 - nu, nu, 0},
    {nu, 1 - nu, 0},
    {0, 0, (1 - 2 * nu) / 2},
  };
  // clang-format on
  solution *= E / ((1 + nu) * (1 - 2 * nu));

  this->computeTangentModuliOnQuad(tangent);
  Real tangent_error = (tangent - solution).norm<L_2>();
  ASSERT_NEAR(tangent_error, 0, 1e-14);

  /* Plane Stress */
  this->plane_stress = true;
  this->updateInternalParameters();
  // clang-format off
  solution = {
    {1, nu, 0},
    {nu, 1, 0},
    {0, 0, (1 - nu) / 2},
  };
  // clang-format on
  solution *= E / (1 - nu * nu);

  this->computeTangentModuliOnQuad(tangent);
  tangent_error = (tangent - solution).norm<L_2>();
  ASSERT_NEAR(tangent_error, 0, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <> void FriendMaterial<MaterialElastic<2>>::testEnergyDensity() {
  Matrix<Real> sigma = {{1, 2}, {2, 4}};
  Matrix<Real> eps = {{1, 0}, {0, 1}};
  Real epot = 0;
  Real solution = 2.5;
  this->computePotentialEnergyOnQuad(eps, sigma, epot);
  ASSERT_NEAR(epot, solution, 1e-14);
}

/* -------------------------------------------------------------------------- */
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
} // namespace
