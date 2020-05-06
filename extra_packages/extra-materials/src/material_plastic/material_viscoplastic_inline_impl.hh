/**
 * @file   material_viscoplastic_inline_impl.hh
 *
 * @author Ramin Aghababaei <ramin.aghababaei@epfl.ch>
 *
 *
 * @brief  Implementation of the inline functions of the material viscoplastic
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

#include "material_viscoplastic.hh"
#include <cmath>

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template <UInt dim>
inline void MaterialViscoPlastic<dim>::computeStressOnQuad(
    const Matrix<Real> & grad_u, const Matrix<Real> & previous_grad_u,
    Matrix<Real> & sigma, const Matrix<Real> & previous_sigma,
    Matrix<Real> & inelastic_strain,
    const Matrix<Real> & previous_inelastic_strain,
    Real & iso_hardening) const {
  // Infinitesimal plastic
  // Compute stress magnitude
  Real s = sigma.doubleDot(sigma);
  Real sigma_mag = sqrt(s);

  // Compute plastic strain increment
  Real factor = (this->ts * this->edot0 * pow(sigma_mag, (this->rate - 1.)) /
                 pow(this->sigma_y + iso_hardening, this->rate));

  Matrix<Real> delta_inelastic_strain(sigma);
  delta_inelastic_strain *= factor;

  // Compute plastic strain increment magnitude
  s = delta_inelastic_strain.doubleDot(delta_inelastic_strain);
  Real dep_mag = std::sqrt(s);

  Matrix<Real> grad_delta_u(grad_u);
  grad_delta_u -= previous_grad_u;

  // Update stress and plastic strain
  Matrix<Real> grad_u_elastic(dim, dim);
  grad_u_elastic = grad_delta_u;
  grad_u_elastic -= delta_inelastic_strain;

  Matrix<Real> sigma_elastic(dim, dim);
  MaterialElastic<dim>::computeStressOnQuad(grad_u_elastic, sigma_elastic);
  sigma += sigma_elastic;

  inelastic_strain += delta_inelastic_strain;

  // Update resistance stress
  iso_hardening = iso_hardening + this->h * dep_mag;

  MaterialPlastic<dim>::computeStressAndInelasticStrainOnQuad(
      grad_delta_u, sigma, previous_sigma, inelastic_strain,
      previous_inelastic_strain, delta_inelastic_strain);
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
inline void MaterialViscoPlastic<dim>::computeTangentModuliOnQuad(
    Matrix<Real> & tangent, const Matrix<Real> & /*grad_u*/,
    const Matrix<Real> & /*previous_grad_u*/,
    const Matrix<Real> & /*sigma_tensor*/,
    const Matrix<Real> & /*previous_sigma_tensor*/,
    const Real & /*iso_hardening*/) const {
  UInt cols = tangent.cols();
  UInt rows = tangent.rows();

  for (UInt m = 0; m < rows; ++m) {
    UInt i = VoigtHelper<dim>::vec[m][0];
    UInt j = VoigtHelper<dim>::vec[m][1];

    for (UInt n = 0; n < cols; ++n) {
      UInt k = VoigtHelper<dim>::vec[n][0];
      UInt l = VoigtHelper<dim>::vec[n][1];
      tangent(m, n) = (i == k) * (j == l) * 2. * this->mu +
                      (i == j) * (k == l) * this->lambda;
      tangent(m, n) -= (m == n) * (m >= dim) * this->mu;
    }
  }
}
