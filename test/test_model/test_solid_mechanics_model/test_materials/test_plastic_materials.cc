/* -------------------------------------------------------------------------- */
#include "material_linear_isotropic_hardening.hh"
#include "solid_mechanics_model.hh"
#include "test_material_fixtures.hh"
#include <gtest/gtest.h>
#include <type_traits>
/* -------------------------------------------------------------------------- */

using namespace akantu;

using types = ::testing::Types<Traits<MaterialLinearIsotropicHardening, 1>,
                               Traits<MaterialLinearIsotropicHardening, 2>,
                               Traits<MaterialLinearIsotropicHardening, 3>>;

/* -------------------------------------------------------------------------- */

template <>
void FriendMaterial<MaterialLinearIsotropicHardening<3>>::testComputeStress() {

  Real E = 1.;
  // Real nu = .3;
  Real nu = 0.;
  Real rho = 1.;
  Real sigma_0 = 1.;
  Real h = 0.;
  Real bulk_modulus_K = E / 3. / (1 - 2. * nu);
  Real shear_modulus_mu = 0.5 * E / (1 + nu);

  setParam("E", E);
  setParam("nu", nu);
  setParam("rho", rho);
  setParam("sigma_y", sigma_0);
  setParam("h", h);

  Matrix<Real> rotation_matrix = getRandomRotation3d();

  Real max_strain = 10.;
  Real strain_steps = 100;
  Real dt = max_strain / strain_steps;
  std::vector<double> steps(strain_steps);
  std::iota(steps.begin(), steps.end(), 0.);

  Matrix<Real> previous_grad_u_rot = Matrix<Real>(3, 3, 0.);
  Matrix<Real> previous_sigma = Matrix<Real>(3, 3, 0.);
  Matrix<Real> previous_sigma_rot = Matrix<Real>(3, 3, 0.);
  Matrix<Real> inelastic_strain_rot = Matrix<Real>(3, 3, 0.);
  Matrix<Real> inelastic_strain = Matrix<Real>(3, 3, 0.);
  Matrix<Real> previous_inelastic_strain = Matrix<Real>(3, 3, 0.);
  Matrix<Real> previous_inelastic_strain_rot = Matrix<Real>(3, 3, 0.);
  Matrix<Real> sigma_rot(3, 3, 0.);
  Real iso_hardening = 0.;
  Real previous_iso_hardening = 0.;

  // hydrostatic loading (should not plastify)
  for (auto && i : steps) {
    auto t = i * dt;

    auto grad_u = this->getHydrostaticStrain(t);
    auto grad_u_rot = this->applyRotation(grad_u, rotation_matrix);

    this->computeStressOnQuad(grad_u_rot, previous_grad_u_rot, sigma_rot,
                              previous_sigma_rot, inelastic_strain_rot,
                              previous_inelastic_strain_rot, iso_hardening,
                              previous_iso_hardening, 0., 0.);

    auto sigma = this->reverseRotation(sigma_rot, rotation_matrix);

    Matrix<Real> sigma_expected =
        t * 3. * bulk_modulus_K * Matrix<Real>::eye(3, 1.);

    Real stress_error = (sigma - sigma_expected).norm<L_inf>();

    ASSERT_NEAR(stress_error, 0., 1e-13);
    ASSERT_NEAR(inelastic_strain_rot.norm<L_inf>(), 0., 1e-13);

    previous_grad_u_rot = grad_u_rot;
    previous_sigma_rot = sigma_rot;
    previous_inelastic_strain_rot = inelastic_strain_rot;
    previous_iso_hardening = iso_hardening;
  }

  // deviatoric loading (should plastify)
  // stress at onset of plastication
  Real beta = sqrt(42);
  Real t_P = sigma_0 / 2. / shear_modulus_mu / beta;
  Matrix<Real> sigma_P = sigma_0 / beta * this->getDeviatoricStrain(1.);

  for (auto && i : steps) {

    auto t = i * dt;
    auto grad_u = this->getDeviatoricStrain(t);
    auto grad_u_rot = this->applyRotation(grad_u, rotation_matrix);
    Real iso_hardening, previous_iso_hardening;

    this->computeStressOnQuad(grad_u_rot, previous_grad_u_rot, sigma_rot,
                              previous_sigma_rot, inelastic_strain_rot,
                              previous_inelastic_strain_rot, iso_hardening,
                              previous_iso_hardening, 0., 0.);

    auto sigma = this->reverseRotation(sigma_rot, rotation_matrix);
    auto inelastic_strain =
        this->reverseRotation(inelastic_strain_rot, rotation_matrix);

    if (t < t_P) {

      Matrix<Real> sigma_expected =
          shear_modulus_mu * (grad_u + grad_u.transpose());

      Real stress_error = (sigma - sigma_expected).norm<L_inf>();
      ASSERT_NEAR(stress_error, 0., 1e-13);
      ASSERT_NEAR(inelastic_strain_rot.norm<L_inf>(), 0., 1e-13);
    } else if (t > t_P + dt) {
      // skip the transition from non plastic to plastic

      auto delta_lambda_expected =
          dt / t * previous_sigma.doubleDot(grad_u + grad_u.transpose()) / 2.;
      auto delta_inelastic_strain_expected =
          delta_lambda_expected * 3. / 2. / sigma_0 * previous_sigma;
      auto inelastic_strain_expected =
          delta_inelastic_strain_expected + previous_inelastic_strain;
      ASSERT_NEAR((inelastic_strain - inelastic_strain_expected).norm<L_inf>(),
                  0., 1e-13);
      auto delta_sigma_expected =
          2. * shear_modulus_mu *
          (0.5 * dt / t * (grad_u + grad_u.transpose()) -
           delta_inelastic_strain_expected);

      auto delta_sigma = sigma - previous_sigma;
      ASSERT_NEAR((delta_sigma_expected - delta_sigma).norm<L_inf>(), 0.,
                  1e-13);
    }
    previous_sigma = sigma;
    previous_sigma_rot = sigma_rot;
    previous_grad_u_rot = grad_u_rot;
    previous_inelastic_strain = inelastic_strain;
    previous_inelastic_strain_rot = inelastic_strain_rot;
  }
}

namespace {

template <typename T>
class TestPlasticMaterialFixture : public ::TestMaterialFixture<T> {};

TYPED_TEST_CASE(TestPlasticMaterialFixture, types);

TYPED_TEST(TestPlasticMaterialFixture, PlasticComputeStress) {
  this->material->testComputeStress();
}
TYPED_TEST(TestPlasticMaterialFixture, PlasticEnergyDensity) {
  this->material->testEnergyDensity();
}
TYPED_TEST(TestPlasticMaterialFixture, PlasticComputeTangentModuli) {
  this->material->testComputeTangentModuli();
}
TYPED_TEST(TestPlasticMaterialFixture, PlasticComputeCelerity) {
  this->material->testCelerity();
}
}

/*****************************************************************/
