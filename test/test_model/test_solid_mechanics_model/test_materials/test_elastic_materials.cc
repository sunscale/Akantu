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

                     Traits<MaterialElasticOrthotropic, 2>,
                     Traits<MaterialElasticOrthotropic, 3>,

                     Traits<MaterialElasticLinearAnisotropic, 2>,
                     Traits<MaterialElasticLinearAnisotropic, 3>>;

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
template <> void FriendMaterial<MaterialElastic<1>>::testCelerity() {
  Real E = 3., rho = 2.;
  setParam("E", E);
  setParam("rho", rho);

  auto wave_speed = this->getCelerity(Element());
  auto solution = std::sqrt(E / rho);
  ASSERT_NEAR(wave_speed, solution, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <> void FriendMaterial<MaterialElastic<2>>::testComputeStress() {
  Real E = 1.;
  Real nu = .3;
  Real sigma_th = 0.3; // thermal stress
  setParam("E", E);
  setParam("nu", nu);

  Real bulk_modulus_K = E / 3. / (1 - 2. * nu);
  Real shear_modulus_mu = 0.5 * E / (1 + nu);

  Matrix<Real> rotation_matrix = getRandomRotation2d();

  auto grad_u = this->getComposedStrain(1.).block(0, 0, 2, 2);

  auto grad_u_rot = this->applyRotation(grad_u, rotation_matrix);

  Matrix<Real> sigma_rot(2, 2);
  this->computeStressOnQuad(grad_u_rot, sigma_rot, sigma_th);

  auto sigma = this->reverseRotation(sigma_rot, rotation_matrix);

  Matrix<Real> identity(2, 2);
  identity.eye();

  Matrix<Real> strain = 0.5 * (grad_u + grad_u.transpose());
  Matrix<Real> deviatoric_strain = strain - 1. / 3. * strain.trace() * identity;

  Matrix<Real> sigma_expected = 2 * shear_modulus_mu * deviatoric_strain +
                                (sigma_th + 2. * bulk_modulus_K) * identity;

  auto diff = sigma - sigma_expected;
  Real stress_error = diff.norm<L_inf>();
  ASSERT_NEAR(stress_error, 0., 1e-14);
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
template <> void FriendMaterial<MaterialElastic<2>>::testCelerity() {
  Real E = 1.;
  Real nu = .3;
  Real rho = 2;
  setParam("E", E);
  setParam("nu", nu);
  setParam("rho", rho);
  auto push_wave_speed = this->getPushWaveSpeed(Element());
  auto celerity = this->getCelerity(Element());

  Real K = E / (3 * (1 - 2 * nu));
  Real mu = E / (2 * (1 + nu));
  Real sol = std::sqrt((K + 4. / 3 * mu) / rho);

  ASSERT_NEAR(push_wave_speed, sol, 1e-14);
  ASSERT_NEAR(celerity, sol, 1e-14);

  auto shear_wave_speed = this->getShearWaveSpeed(Element());

  sol = std::sqrt(mu / rho);

  ASSERT_NEAR(shear_wave_speed, sol, 1e-14);
}

/* -------------------------------------------------------------------------- */

template <> void FriendMaterial<MaterialElastic<3>>::testComputeStress() {
  Real E = 1.;
  Real nu = .3;
  Real sigma_th = 0.3; // thermal stress
  setParam("E", E);
  setParam("nu", nu);

  Real bulk_modulus_K = E / 3. / (1 - 2. * nu);
  Real shear_modulus_mu = 0.5 * E / (1 + nu);

  Matrix<Real> rotation_matrix = getRandomRotation3d();

  auto grad_u = this->getComposedStrain(1.);

  auto grad_u_rot = this->applyRotation(grad_u, rotation_matrix);

  Matrix<Real> sigma_rot(3, 3);
  this->computeStressOnQuad(grad_u_rot, sigma_rot, sigma_th);

  auto sigma = this->reverseRotation(sigma_rot, rotation_matrix);

  Matrix<Real> identity(3, 3);
  identity.eye();

  Matrix<Real> strain = 0.5 * (grad_u + grad_u.transpose());
  Matrix<Real> deviatoric_strain = strain - 1. / 3. * strain.trace() * identity;

  Matrix<Real> sigma_expected = 2 * shear_modulus_mu * deviatoric_strain +
                                (sigma_th + 3. * bulk_modulus_K) * identity;

  auto diff = sigma - sigma_expected;
  Real stress_error = diff.norm<L_inf>();
  ASSERT_NEAR(stress_error, 0., 1e-14);
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
template <> void FriendMaterial<MaterialElastic<3>>::testCelerity() {
  Real E = 1.;
  Real nu = .3;
  Real rho = 2;
  setParam("E", E);
  setParam("nu", nu);
  setParam("rho", rho);

  auto push_wave_speed = this->getPushWaveSpeed(Element());
  auto celerity = this->getCelerity(Element());

  Real K = E / (3 * (1 - 2 * nu));
  Real mu = E / (2 * (1 + nu));
  Real sol = std::sqrt((K + 4. / 3 * mu) / rho);

  ASSERT_NEAR(push_wave_speed, sol, 1e-14);
  ASSERT_NEAR(celerity, sol, 1e-14);

  auto shear_wave_speed = this->getShearWaveSpeed(Element());

  sol = std::sqrt(mu / rho);

  ASSERT_NEAR(shear_wave_speed, sol, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<MaterialElasticOrthotropic<2>>::testComputeStress() {
  UInt Dim = 2;

  this->E1 = 1.;
  this->E2 = 2.;
  this->nu12 = 0.1;
  this->G12 = 2.;

  // material frame of reference is rotate by rotation_matrix starting from
  // canonical basis
  Matrix<Real> rotation_matrix = getRandomRotation2d();

  // canonical basis as expressed in the material frame of reference, as
  // required by MaterialElasticOrthotropic class (it is simply given by the
  // columns of the rotation_matrix; the lines give the material basis expressed
  // in the canonical frame of reference)
  *this->dir_vecs[0] = rotation_matrix(0);
  *this->dir_vecs[1] = rotation_matrix(1);

  // set internal Cijkl matrix expressed in the canonical frame of reference
  this->updateInternalParameters();

  // gradient in material frame of reference
  auto grad_u = this->getComposedStrain(2.).block(0, 0, 2, 2);

  // gradient in canonical basis (we need to rotate *back* to the canonical
  // basis)
  auto grad_u_rot = this->reverseRotation(grad_u, rotation_matrix);

  // stress in the canonical basis
  Matrix<Real> sigma_rot(2, 2);
  this->computeStressOnQuad(grad_u_rot, sigma_rot);

  // stress in the material reference (we need to apply the rotation)
  auto sigma = this->applyRotation(sigma_rot, rotation_matrix);

  // construction of Cijkl engineering tensor in the *material* frame of
  // reference
  // ref: http://solidmechanics.org/Text/Chapter3_2/Chapter3_2.php#Sect3_2_13
  Real nu21 = nu12 * E2 / E1;
  Real gamma = 1 / (1 - nu12 * nu21);

  Matrix<Real> C_expected(2 * Dim, 2 * Dim, 0);
  C_expected(0, 0) = gamma * E1;
  C_expected(1, 1) = gamma * E2;
  C_expected(2, 2) = G12;

  C_expected(1, 0) = C_expected(0, 1) = gamma * E1 * nu21;

  // epsilon is computed directly in the *material* frame of reference
  Matrix<Real> epsilon = 0.5 * (grad_u + grad_u.transpose());

  // sigma_expected is computed directly in the *material* frame of reference
  Matrix<Real> sigma_expected(Dim, Dim);
  for (UInt i = 0; i < Dim; ++i) {
    for (UInt j = 0; j < Dim; ++j) {
      sigma_expected(i, i) += C_expected(i, j) * epsilon(j, j);
    }
  }

  sigma_expected(0, 1) = sigma_expected(1, 0) =
      C_expected(2, 2) * 2 * epsilon(0, 1);

  // sigmas are checked in the *material* frame of reference
  auto diff = sigma - sigma_expected;
  Real stress_error = diff.norm<L_inf>();
  ASSERT_NEAR(stress_error, 0., 1e-13);
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<MaterialElasticOrthotropic<2>>::testEnergyDensity() {
  Matrix<Real> sigma = {{1, 2}, {2, 4}};
  Matrix<Real> eps = {{1, 0}, {0, 1}};
  Real epot = 0;
  Real solution = 2.5;
  this->computePotentialEnergyOnQuad(eps, sigma, epot);
  ASSERT_NEAR(epot, solution, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<MaterialElasticOrthotropic<2>>::testComputeTangentModuli() {

  // Note: for this test material and canonical basis coincide
  Vector<Real> n1 = {1, 0};
  Vector<Real> n2 = {0, 1};
  Real E1 = 1.;
  Real E2 = 2.;
  Real nu12 = 0.1;
  Real G12 = 2.;

  *this->dir_vecs[0] = n1;
  *this->dir_vecs[1] = n2;

  this->E1 = E1;
  this->E2 = E2;
  this->nu12 = nu12;
  this->G12 = G12;

  // set internal Cijkl matrix expressed in the canonical frame of reference
  this->updateInternalParameters();

  // construction of Cijkl engineering tensor in the *material* frame of
  // reference
  // ref: http://solidmechanics.org/Text/Chapter3_2/Chapter3_2.php#Sect3_2_13
  Real nu21 = nu12 * E2 / E1;
  Real gamma = 1 / (1 - nu12 * nu21);

  Matrix<Real> C_expected(3, 3);
  C_expected(0, 0) = gamma * E1;
  C_expected(1, 1) = gamma * E2;
  C_expected(2, 2) = G12;

  C_expected(1, 0) = C_expected(0, 1) = gamma * E1 * nu21;

  Matrix<Real> tangent(3, 3);
  this->computeTangentModuliOnQuad(tangent);

  Real tangent_error = (tangent - C_expected).norm<L_2>();
  ASSERT_NEAR(tangent_error, 0, 1e-14);
}

/* -------------------------------------------------------------------------- */

template <> void FriendMaterial<MaterialElasticOrthotropic<2>>::testCelerity() {

  // Note: for this test material and canonical basis coincide
  Vector<Real> n1 = {1, 0};
  Vector<Real> n2 = {0, 1};
  Real E1 = 1.;
  Real E2 = 2.;
  Real nu12 = 0.1;
  Real G12 = 2.;
  Real rho = 2.5;

  setParamNoUpdate("n1", n1);
  setParamNoUpdate("n2", n2);
  setParamNoUpdate("E1", E1);
  setParamNoUpdate("E2", E2);
  setParamNoUpdate("nu12", nu12);
  setParamNoUpdate("G12", G12);
  setParamNoUpdate("rho", rho);

  // set internal Cijkl matrix expressed in the canonical frame of reference
  this->updateInternalParameters();

  // construction of Cijkl engineering tensor in the *material* frame of
  // reference
  // ref: http://solidmechanics.org/Text/Chapter3_2/Chapter3_2.php#Sect3_2_13
  Real nu21 = nu12 * E2 / E1;
  Real gamma = 1 / (1 - nu12 * nu21);

  Matrix<Real> C_expected(3, 3);
  C_expected(0, 0) = gamma * E1;
  C_expected(1, 1) = gamma * E2;
  C_expected(2, 2) = G12;

  C_expected(1, 0) = C_expected(0, 1) = gamma * E1 * nu21;

  Vector<Real> eig_expected(3);
  C_expected.eig(eig_expected);

  auto celerity_expected = std::sqrt(eig_expected(0) / rho);

  auto celerity = this->getCelerity(Element());

  ASSERT_NEAR(celerity_expected, celerity, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<MaterialElasticOrthotropic<3>>::testComputeStress() {
  UInt Dim = 3;
  Real E1 = 1.;
  Real E2 = 2.;
  Real E3 = 3.;
  Real nu12 = 0.1;
  Real nu13 = 0.2;
  Real nu23 = 0.3;
  Real G12 = 2.;
  Real G13 = 3.;
  Real G23 = 1.;

  this->E1 = E1;
  this->E2 = E2;
  this->E3 = E3;
  this->nu12 = nu12;
  this->nu13 = nu13;
  this->nu23 = nu23;
  this->G12 = G12;
  this->G13 = G13;
  this->G23 = G23;

  // material frame of reference is rotate by rotation_matrix starting from
  // canonical basis
  Matrix<Real> rotation_matrix = getRandomRotation3d();

  // canonical basis as expressed in the material frame of reference, as
  // required by MaterialElasticOrthotropic class (it is simply given by the
  // columns of the rotation_matrix; the lines give the material basis expressed
  // in the canonical frame of reference)
  *this->dir_vecs[0] = rotation_matrix(0);
  *this->dir_vecs[1] = rotation_matrix(1);
  *this->dir_vecs[2] = rotation_matrix(2);

  // set internal Cijkl matrix expressed in the canonical frame of reference
  this->updateInternalParameters();

  // gradient in material frame of reference
  auto grad_u = this->getComposedStrain(2.);

  // gradient in canonical basis (we need to rotate *back* to the canonical
  // basis)
  auto grad_u_rot = this->reverseRotation(grad_u, rotation_matrix);

  // stress in the canonical basis
  Matrix<Real> sigma_rot(3, 3);
  this->computeStressOnQuad(grad_u_rot, sigma_rot);

  // stress in the material reference (we need to apply the rotation)
  auto sigma = this->applyRotation(sigma_rot, rotation_matrix);

  // construction of Cijkl engineering tensor in the *material* frame of
  // reference
  // ref: http://solidmechanics.org/Text/Chapter3_2/Chapter3_2.php#Sect3_2_13
  Real nu21 = nu12 * E2 / E1;
  Real nu31 = nu13 * E3 / E1;
  Real nu32 = nu23 * E3 / E2;
  Real gamma = 1 / (1 - nu12 * nu21 - nu23 * nu32 - nu31 * nu13 -
                    2 * nu21 * nu32 * nu13);

  Matrix<Real> C_expected(6, 6);
  C_expected(0, 0) = gamma * E1 * (1 - nu23 * nu32);
  C_expected(1, 1) = gamma * E2 * (1 - nu13 * nu31);
  C_expected(2, 2) = gamma * E3 * (1 - nu12 * nu21);

  C_expected(1, 0) = C_expected(0, 1) = gamma * E1 * (nu21 + nu31 * nu23);
  C_expected(2, 0) = C_expected(0, 2) = gamma * E1 * (nu31 + nu21 * nu32);
  C_expected(2, 1) = C_expected(1, 2) = gamma * E2 * (nu32 + nu12 * nu31);

  C_expected(3, 3) = G23;
  C_expected(4, 4) = G13;
  C_expected(5, 5) = G12;

  // epsilon is computed directly in the *material* frame of reference
  Matrix<Real> epsilon = 0.5 * (grad_u + grad_u.transpose());

  // sigma_expected is computed directly in the *material* frame of reference
  Matrix<Real> sigma_expected(Dim, Dim);
  for (UInt i = 0; i < Dim; ++i) {
    for (UInt j = 0; j < Dim; ++j) {
      sigma_expected(i, i) += C_expected(i, j) * epsilon(j, j);
    }
  }

  sigma_expected(0, 1) = C_expected(5, 5) * 2 * epsilon(0, 1);
  sigma_expected(0, 2) = C_expected(4, 4) * 2 * epsilon(0, 2);
  sigma_expected(1, 2) = C_expected(3, 3) * 2 * epsilon(1, 2);
  sigma_expected(1, 0) = sigma_expected(0, 1);
  sigma_expected(2, 0) = sigma_expected(0, 2);
  sigma_expected(2, 1) = sigma_expected(1, 2);

  // sigmas are checked in the *material* frame of reference
  auto diff = sigma - sigma_expected;
  Real stress_error = diff.norm<L_inf>();
  ASSERT_NEAR(stress_error, 0., 1e-13);
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<MaterialElasticOrthotropic<3>>::testEnergyDensity() {
  Matrix<Real> sigma = {{1, 2, 3}, {2, 4, 5}, {3, 5, 6}};
  Matrix<Real> eps = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  Real epot = 0;
  Real solution = 5.5;
  this->computePotentialEnergyOnQuad(eps, sigma, epot);
  ASSERT_NEAR(epot, solution, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<MaterialElasticOrthotropic<3>>::testComputeTangentModuli() {

  // Note: for this test material and canonical basis coincide
  UInt Dim = 3;
  Vector<Real> n1 = {1, 0, 0};
  Vector<Real> n2 = {0, 1, 0};
  Vector<Real> n3 = {0, 0, 1};
  Real E1 = 1.;
  Real E2 = 2.;
  Real E3 = 3.;
  Real nu12 = 0.1;
  Real nu13 = 0.2;
  Real nu23 = 0.3;
  Real G12 = 2.;
  Real G13 = 3.;
  Real G23 = 1.;

  *this->dir_vecs[0] = n1;
  *this->dir_vecs[1] = n2;
  *this->dir_vecs[2] = n3;
  this->E1 = E1;
  this->E2 = E2;
  this->E3 = E3;
  this->nu12 = nu12;
  this->nu13 = nu13;
  this->nu23 = nu23;
  this->G12 = G12;
  this->G13 = G13;
  this->G23 = G23;

  // set internal Cijkl matrix expressed in the canonical frame of reference
  this->updateInternalParameters();

  // construction of Cijkl engineering tensor in the *material* frame of
  // reference
  // ref: http://solidmechanics.org/Text/Chapter3_2/Chapter3_2.php#Sect3_2_13
  Real nu21 = nu12 * E2 / E1;
  Real nu31 = nu13 * E3 / E1;
  Real nu32 = nu23 * E3 / E2;
  Real gamma = 1 / (1 - nu12 * nu21 - nu23 * nu32 - nu31 * nu13 -
                    2 * nu21 * nu32 * nu13);

  Matrix<Real> C_expected(2 * Dim, 2 * Dim, 0);
  C_expected(0, 0) = gamma * E1 * (1 - nu23 * nu32);
  C_expected(1, 1) = gamma * E2 * (1 - nu13 * nu31);
  C_expected(2, 2) = gamma * E3 * (1 - nu12 * nu21);

  C_expected(1, 0) = C_expected(0, 1) = gamma * E1 * (nu21 + nu31 * nu23);
  C_expected(2, 0) = C_expected(0, 2) = gamma * E1 * (nu31 + nu21 * nu32);
  C_expected(2, 1) = C_expected(1, 2) = gamma * E2 * (nu32 + nu12 * nu31);

  C_expected(3, 3) = G23;
  C_expected(4, 4) = G13;
  C_expected(5, 5) = G12;

  Matrix<Real> tangent(6, 6);
  this->computeTangentModuliOnQuad(tangent);

  Real tangent_error = (tangent - C_expected).norm<L_2>();
  ASSERT_NEAR(tangent_error, 0, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <> void FriendMaterial<MaterialElasticOrthotropic<3>>::testCelerity() {

  // Note: for this test material and canonical basis coincide
  UInt Dim = 3;
  Vector<Real> n1 = {1, 0, 0};
  Vector<Real> n2 = {0, 1, 0};
  Vector<Real> n3 = {0, 0, 1};
  Real E1 = 1.;
  Real E2 = 2.;
  Real E3 = 3.;
  Real nu12 = 0.1;
  Real nu13 = 0.2;
  Real nu23 = 0.3;
  Real G12 = 2.;
  Real G13 = 3.;
  Real G23 = 1.;
  Real rho = 2.3;

  setParamNoUpdate("n1", n1);
  setParamNoUpdate("n2", n2);
  setParamNoUpdate("n3", n3);
  setParamNoUpdate("E1", E1);
  setParamNoUpdate("E2", E2);
  setParamNoUpdate("E3", E3);
  setParamNoUpdate("nu12", nu12);
  setParamNoUpdate("nu13", nu13);
  setParamNoUpdate("nu23", nu23);
  setParamNoUpdate("G12", G12);
  setParamNoUpdate("G13", G13);
  setParamNoUpdate("G23", G23);
  setParamNoUpdate("rho", rho);

  // set internal Cijkl matrix expressed in the canonical frame of reference
  this->updateInternalParameters();

  // construction of Cijkl engineering tensor in the *material* frame of
  // reference
  // ref: http://solidmechanics.org/Text/Chapter3_2/Chapter3_2.php#Sect3_2_13
  Real nu21 = nu12 * E2 / E1;
  Real nu31 = nu13 * E3 / E1;
  Real nu32 = nu23 * E3 / E2;
  Real gamma = 1 / (1 - nu12 * nu21 - nu23 * nu32 - nu31 * nu13 -
                    2 * nu21 * nu32 * nu13);

  Matrix<Real> C_expected(2 * Dim, 2 * Dim, 0);
  C_expected(0, 0) = gamma * E1 * (1 - nu23 * nu32);
  C_expected(1, 1) = gamma * E2 * (1 - nu13 * nu31);
  C_expected(2, 2) = gamma * E3 * (1 - nu12 * nu21);

  C_expected(1, 0) = C_expected(0, 1) = gamma * E1 * (nu21 + nu31 * nu23);
  C_expected(2, 0) = C_expected(0, 2) = gamma * E1 * (nu31 + nu21 * nu32);
  C_expected(2, 1) = C_expected(1, 2) = gamma * E2 * (nu32 + nu12 * nu31);

  C_expected(3, 3) = G23;
  C_expected(4, 4) = G13;
  C_expected(5, 5) = G12;

  Vector<Real> eig_expected(6);
  C_expected.eig(eig_expected);

  auto celerity_expected = std::sqrt(eig_expected(0) / rho);

  auto celerity = this->getCelerity(Element());

  ASSERT_NEAR(celerity_expected, celerity, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<MaterialElasticLinearAnisotropic<2>>::testComputeStress() {
  UInt Dim = 2;
  Matrix<Real> C = {
      {1.0, 0.3, 0.4}, {0.3, 2.0, 0.1}, {0.4, 0.1, 1.5},
  };

  setParamNoUpdate("C11", C(0, 0));
  setParamNoUpdate("C12", C(0, 1));
  setParamNoUpdate("C13", C(0, 2));
  setParamNoUpdate("C22", C(1, 1));
  setParamNoUpdate("C23", C(1, 2));
  setParamNoUpdate("C33", C(2, 2));

  // material frame of reference is rotate by rotation_matrix starting from
  // canonical basis
  Matrix<Real> rotation_matrix = getRandomRotation2d();

  // canonical basis as expressed in the material frame of reference, as
  // required by MaterialElasticLinearAnisotropic class (it is simply given by
  // the columns of the rotation_matrix; the lines give the material basis
  // expressed in the canonical frame of reference)
  *this->dir_vecs[0] = rotation_matrix(0);
  *this->dir_vecs[1] = rotation_matrix(1);

  // set internal Cijkl matrix expressed in the canonical frame of reference
  this->updateInternalParameters();

  // gradient in material frame of reference
  auto grad_u = this->getComposedStrain(1.).block(0, 0, 2, 2);

  // gradient in canonical basis (we need to rotate *back* to the canonical
  // basis)
  auto grad_u_rot = this->reverseRotation(grad_u, rotation_matrix);

  // stress in the canonical basis
  Matrix<Real> sigma_rot(2, 2);
  this->computeStressOnQuad(grad_u_rot, sigma_rot);

  // stress in the material reference (we need to apply the rotation)
  auto sigma = this->applyRotation(sigma_rot, rotation_matrix);

  // epsilon is computed directly in the *material* frame of reference
  Matrix<Real> epsilon = 0.5 * (grad_u + grad_u.transpose());
  Vector<Real> epsilon_voigt(3);
  epsilon_voigt(0) = epsilon(0, 0);
  epsilon_voigt(1) = epsilon(1, 1);
  epsilon_voigt(2) = 2 * epsilon(0, 1);

  // sigma_expected is computed directly in the *material* frame of reference
  Vector<Real> sigma_voigt = C * epsilon_voigt;
  Matrix<Real> sigma_expected(Dim, Dim);
  sigma_expected(0, 0) = sigma_voigt(0);
  sigma_expected(1, 1) = sigma_voigt(1);
  sigma_expected(0, 1) = sigma_expected(1, 0) = sigma_voigt(2);

  // sigmas are checked in the *material* frame of reference
  auto diff = sigma - sigma_expected;
  Real stress_error = diff.norm<L_inf>();
  ASSERT_NEAR(stress_error, 0., 1e-13);
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<MaterialElasticLinearAnisotropic<2>>::testEnergyDensity() {
  Matrix<Real> sigma = {{1, 2}, {2, 4}};
  Matrix<Real> eps = {{1, 0}, {0, 1}};
  Real epot = 0;
  Real solution = 2.5;
  this->computePotentialEnergyOnQuad(eps, sigma, epot);
  ASSERT_NEAR(epot, solution, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<
    MaterialElasticLinearAnisotropic<2>>::testComputeTangentModuli() {

  // Note: for this test material and canonical basis coincide
  Matrix<Real> C = {
      {1.0, 0.3, 0.4}, {0.3, 2.0, 0.1}, {0.4, 0.1, 1.5},
  };

  setParamNoUpdate("C11", C(0, 0));
  setParamNoUpdate("C12", C(0, 1));
  setParamNoUpdate("C13", C(0, 2));
  setParamNoUpdate("C22", C(1, 1));
  setParamNoUpdate("C23", C(1, 2));
  setParamNoUpdate("C33", C(2, 2));

  // set internal Cijkl matrix expressed in the canonical frame of reference
  this->updateInternalParameters();

  Matrix<Real> tangent(3, 3);
  this->computeTangentModuliOnQuad(tangent);

  Real tangent_error = (tangent - C).norm<L_2>();
  ASSERT_NEAR(tangent_error, 0, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<MaterialElasticLinearAnisotropic<2>>::testCelerity() {

  // Note: for this test material and canonical basis coincide
  Matrix<Real> C = {
      {1.0, 0.3, 0.4}, {0.3, 2.0, 0.1}, {0.4, 0.1, 1.5},
  };

  Real rho = 2.7;

  setParamNoUpdate("C11", C(0, 0));
  setParamNoUpdate("C12", C(0, 1));
  setParamNoUpdate("C13", C(0, 2));
  setParamNoUpdate("C22", C(1, 1));
  setParamNoUpdate("C23", C(1, 2));
  setParamNoUpdate("C33", C(2, 2));
  setParamNoUpdate("rho", rho);

  // set internal Cijkl matrix expressed in the canonical frame of reference
  this->updateInternalParameters();

  Vector<Real> eig_expected(3);
  C.eig(eig_expected);

  auto celerity_expected = std::sqrt(eig_expected(0) / rho);

  auto celerity = this->getCelerity(Element());

  ASSERT_NEAR(celerity_expected, celerity, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<MaterialElasticLinearAnisotropic<3>>::testComputeStress() {
  UInt Dim = 3;
  Matrix<Real> C = {
      {1.0, 0.3, 0.4, 0.3, 0.2, 0.1}, {0.3, 2.0, 0.1, 0.2, 0.3, 0.2},
      {0.4, 0.1, 1.5, 0.1, 0.4, 0.3}, {0.3, 0.2, 0.1, 2.4, 0.1, 0.4},
      {0.2, 0.3, 0.4, 0.1, 0.9, 0.1}, {0.1, 0.2, 0.3, 0.4, 0.1, 1.2},
  };

  setParamNoUpdate("C11", C(0, 0));
  setParamNoUpdate("C12", C(0, 1));
  setParamNoUpdate("C13", C(0, 2));
  setParamNoUpdate("C14", C(0, 3));
  setParamNoUpdate("C15", C(0, 4));
  setParamNoUpdate("C16", C(0, 5));
  setParamNoUpdate("C22", C(1, 1));
  setParamNoUpdate("C23", C(1, 2));
  setParamNoUpdate("C24", C(1, 3));
  setParamNoUpdate("C25", C(1, 4));
  setParamNoUpdate("C26", C(1, 5));
  setParamNoUpdate("C33", C(2, 2));
  setParamNoUpdate("C34", C(2, 3));
  setParamNoUpdate("C35", C(2, 4));
  setParamNoUpdate("C36", C(2, 5));
  setParamNoUpdate("C44", C(3, 3));
  setParamNoUpdate("C45", C(3, 4));
  setParamNoUpdate("C46", C(3, 5));
  setParamNoUpdate("C55", C(4, 4));
  setParamNoUpdate("C56", C(4, 5));
  setParamNoUpdate("C66", C(5, 5));

  // material frame of reference is rotate by rotation_matrix starting from
  // canonical basis
  Matrix<Real> rotation_matrix = getRandomRotation3d();

  // canonical basis as expressed in the material frame of reference, as
  // required by MaterialElasticLinearAnisotropic class (it is simply given by
  // the columns of the rotation_matrix; the lines give the material basis
  // expressed in the canonical frame of reference)
  *this->dir_vecs[0] = rotation_matrix(0);
  *this->dir_vecs[1] = rotation_matrix(1);
  *this->dir_vecs[2] = rotation_matrix(2);

  // set internal Cijkl matrix expressed in the canonical frame of reference
  this->updateInternalParameters();

  // gradient in material frame of reference
  auto grad_u = this->getComposedStrain(2.);

  // gradient in canonical basis (we need to rotate *back* to the canonical
  // basis)
  auto grad_u_rot = this->reverseRotation(grad_u, rotation_matrix);

  // stress in the canonical basis
  Matrix<Real> sigma_rot(3, 3);
  this->computeStressOnQuad(grad_u_rot, sigma_rot);

  // stress in the material reference (we need to apply the rotation)
  auto sigma = this->applyRotation(sigma_rot, rotation_matrix);

  // epsilon is computed directly in the *material* frame of reference
  Matrix<Real> epsilon = 0.5 * (grad_u + grad_u.transpose());
  Vector<Real> epsilon_voigt(6);
  epsilon_voigt(0) = epsilon(0, 0);
  epsilon_voigt(1) = epsilon(1, 1);
  epsilon_voigt(2) = epsilon(2, 2);
  epsilon_voigt(3) = 2 * epsilon(1, 2);
  epsilon_voigt(4) = 2 * epsilon(0, 2);
  epsilon_voigt(5) = 2 * epsilon(0, 1);

  // sigma_expected is computed directly in the *material* frame of reference
  Vector<Real> sigma_voigt = C * epsilon_voigt;
  Matrix<Real> sigma_expected(Dim, Dim);
  sigma_expected(0, 0) = sigma_voigt(0);
  sigma_expected(1, 1) = sigma_voigt(1);
  sigma_expected(2, 2) = sigma_voigt(2);
  sigma_expected(1, 2) = sigma_expected(2, 1) = sigma_voigt(3);
  sigma_expected(0, 2) = sigma_expected(2, 0) = sigma_voigt(4);
  sigma_expected(0, 1) = sigma_expected(1, 0) = sigma_voigt(5);

  // sigmas are checked in the *material* frame of reference
  auto diff = sigma - sigma_expected;
  Real stress_error = diff.norm<L_inf>();
  ASSERT_NEAR(stress_error, 0., 1e-13);
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<MaterialElasticLinearAnisotropic<3>>::testEnergyDensity() {
  Matrix<Real> sigma = {{1, 2, 3}, {2, 4, 5}, {3, 5, 6}};
  Matrix<Real> eps = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  Real epot = 0;
  Real solution = 5.5;
  this->computePotentialEnergyOnQuad(eps, sigma, epot);
  ASSERT_NEAR(epot, solution, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<
    MaterialElasticLinearAnisotropic<3>>::testComputeTangentModuli() {

  // Note: for this test material and canonical basis coincide
  Matrix<Real> C = {
      {1.0, 0.3, 0.4, 0.3, 0.2, 0.1}, {0.3, 2.0, 0.1, 0.2, 0.3, 0.2},
      {0.4, 0.1, 1.5, 0.1, 0.4, 0.3}, {0.3, 0.2, 0.1, 2.4, 0.1, 0.4},
      {0.2, 0.3, 0.4, 0.1, 0.9, 0.1}, {0.1, 0.2, 0.3, 0.4, 0.1, 1.2},
  };

  setParamNoUpdate("C11", C(0, 0));
  setParamNoUpdate("C12", C(0, 1));
  setParamNoUpdate("C13", C(0, 2));
  setParamNoUpdate("C14", C(0, 3));
  setParamNoUpdate("C15", C(0, 4));
  setParamNoUpdate("C16", C(0, 5));
  setParamNoUpdate("C22", C(1, 1));
  setParamNoUpdate("C23", C(1, 2));
  setParamNoUpdate("C24", C(1, 3));
  setParamNoUpdate("C25", C(1, 4));
  setParamNoUpdate("C26", C(1, 5));
  setParamNoUpdate("C33", C(2, 2));
  setParamNoUpdate("C34", C(2, 3));
  setParamNoUpdate("C35", C(2, 4));
  setParamNoUpdate("C36", C(2, 5));
  setParamNoUpdate("C44", C(3, 3));
  setParamNoUpdate("C45", C(3, 4));
  setParamNoUpdate("C46", C(3, 5));
  setParamNoUpdate("C55", C(4, 4));
  setParamNoUpdate("C56", C(4, 5));
  setParamNoUpdate("C66", C(5, 5));

  // set internal Cijkl matrix expressed in the canonical frame of reference
  this->updateInternalParameters();

  Matrix<Real> tangent(6, 6);
  this->computeTangentModuliOnQuad(tangent);

  Real tangent_error = (tangent - C).norm<L_2>();
  ASSERT_NEAR(tangent_error, 0, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<MaterialElasticLinearAnisotropic<3>>::testCelerity() {

  // Note: for this test material and canonical basis coincide
  Matrix<Real> C = {
      {1.0, 0.3, 0.4, 0.3, 0.2, 0.1}, {0.3, 2.0, 0.1, 0.2, 0.3, 0.2},
      {0.4, 0.1, 1.5, 0.1, 0.4, 0.3}, {0.3, 0.2, 0.1, 2.4, 0.1, 0.4},
      {0.2, 0.3, 0.4, 0.1, 0.9, 0.1}, {0.1, 0.2, 0.3, 0.4, 0.1, 1.2},
  };
  Real rho = 2.9;

  setParamNoUpdate("C11", C(0, 0));
  setParamNoUpdate("C12", C(0, 1));
  setParamNoUpdate("C13", C(0, 2));
  setParamNoUpdate("C14", C(0, 3));
  setParamNoUpdate("C15", C(0, 4));
  setParamNoUpdate("C16", C(0, 5));
  setParamNoUpdate("C22", C(1, 1));
  setParamNoUpdate("C23", C(1, 2));
  setParamNoUpdate("C24", C(1, 3));
  setParamNoUpdate("C25", C(1, 4));
  setParamNoUpdate("C26", C(1, 5));
  setParamNoUpdate("C33", C(2, 2));
  setParamNoUpdate("C34", C(2, 3));
  setParamNoUpdate("C35", C(2, 4));
  setParamNoUpdate("C36", C(2, 5));
  setParamNoUpdate("C44", C(3, 3));
  setParamNoUpdate("C45", C(3, 4));
  setParamNoUpdate("C46", C(3, 5));
  setParamNoUpdate("C55", C(4, 4));
  setParamNoUpdate("C56", C(4, 5));
  setParamNoUpdate("C66", C(5, 5));
  setParamNoUpdate("rho", rho);

  // set internal Cijkl matrix expressed in the canonical frame of reference
  this->updateInternalParameters();

  Vector<Real> eig_expected(6);
  C.eig(eig_expected);

  auto celerity_expected = std::sqrt(eig_expected(0) / rho);

  auto celerity = this->getCelerity(Element());

  ASSERT_NEAR(celerity_expected, celerity, 1e-14);
}

/* -------------------------------------------------------------------------- */

namespace {

template <typename T>
class TestElasticMaterialFixture : public ::TestMaterialFixture<T> {};

TYPED_TEST_CASE(TestElasticMaterialFixture, types);

TYPED_TEST(TestElasticMaterialFixture, ComputeStress) {
  this->material->testComputeStress();
}

TYPED_TEST(TestElasticMaterialFixture, EnergyDensity) {
  this->material->testEnergyDensity();
}

TYPED_TEST(TestElasticMaterialFixture, ComputeTangentModuli) {
  this->material->testComputeTangentModuli();
}

TYPED_TEST(TestElasticMaterialFixture, ElasticComputeCelerity) {
  this->material->testCelerity();
}

} // namespace
