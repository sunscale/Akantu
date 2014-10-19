/**
 * @file   material_linear_isotropic_hardening_inline_impl.cc
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Ramin Aghababaei <ramin.aghababaei@epfl.ch>
 *
 * @date creation: Thu Oct 03 2013
 * @date last modification: Mon Jun 09 2014
 *
 * @brief  Implementation of the inline functions of the material plasticity
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "material_linear_isotropic_hardening.hh"


/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template<UInt dim>
inline void
MaterialLinearIsotropicHardening<dim>::computeStressOnQuad(const Matrix<Real> & grad_u,
                                                           const Matrix<Real> & previous_grad_u,
                                                           Matrix<Real> & sigma,
                                                           const Matrix<Real> & previous_sigma,
                                                           Matrix<Real> & inelastic_strain,
                                                           const Matrix<Real> & previous_inelastic_strain,
                                                           Real & iso_hardening,
                                                           const Real & previous_iso_hardening,
                                                           const Real & sigma_th,
                                                           const Real & previous_sigma_th) {
  //Infinitesimal plasticity
  //Real r=iso_hardening;
  Real dp=0.0;
  Real d_dp=0.0;
  UInt n=0;

  Real delta_sigma_th = sigma_th - previous_sigma_th;

  Matrix<Real> grad_delta_u(grad_u);
  grad_delta_u -= previous_grad_u;

  //Compute trial stress, sigma_tr
  Matrix<Real> sigma_tr(dim, dim);
  MaterialElastic<dim>::computeStressOnQuad(grad_delta_u, sigma_tr, delta_sigma_th);
  sigma_tr += previous_sigma;

  // Compute deviatoric trial stress,  sigma_tr_dev
  Matrix<Real> sigma_tr_dev(sigma_tr);
  sigma_tr_dev -= Matrix<Real>::eye(dim, sigma_tr.trace() / 3.0);

  // Compute effective deviatoric trial stress
  Real s = sigma_tr_dev.doubleDot(sigma_tr_dev);
  Real sigma_tr_dev_eff = std::sqrt(3./2. * s);

  const Real iso_hardening_t = previous_iso_hardening;
  //Loop for correcting stress based on yield function
  while ((sigma_tr_dev_eff - iso_hardening - this->sigma_y) > 0) {

    d_dp = (sigma_tr_dev_eff - 3. * this->mu *dp -  iso_hardening - this->sigma_y)
      / (3. * this->mu + this->h);

    //r = r +  h * dp;
    iso_hardening = iso_hardening_t + this->h * d_dp;
    dp = dp + d_dp;

    ++n;

    /// TODO : explicit this criterion with an error message
    if ((d_dp < 1e-5) || (n>50))
      break;
  }

  //Update internal variable
  Matrix<Real> delta_inelastic_strain(dim, dim, 0.);
  if (std::abs(sigma_tr_dev_eff) >
      sigma_tr_dev.norm<L_inf>() * Math::getTolerance()) {
    delta_inelastic_strain.copy(sigma_tr_dev);
    delta_inelastic_strain *= 3./2. * dp / sigma_tr_dev_eff;
  }

  MaterialPlastic<dim>::computeStressAndInelasticStrainOnQuad(grad_delta_u, sigma, previous_sigma,
                                                              inelastic_strain, previous_inelastic_strain,
                                                              delta_inelastic_strain);
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
inline void
MaterialLinearIsotropicHardening<dim>::computeTangentModuliOnQuad(Matrix<Real> & tangent,
                                                                  const Matrix<Real> & grad_u,
                                                                  const Matrix<Real> & previous_grad_u,
                                                                  const Matrix<Real> & sigma_tensor,
                                                                  const Matrix<Real> & previous_sigma_tensor,
                                                                  const Real & iso_hardening) const {
  // Real r=iso_hardening;

  // Matrix<Real> grad_delta_u(grad_u);
  // grad_delta_u -= previous_grad_u;

  // //Compute trial stress, sigma_tr
  // Matrix<Real> sigma_tr(dim, dim);
  // MaterialElastic<dim>::computeStressOnQuad(grad_delta_u, sigma_tr);
  // sigma_tr += previous_sigma_tensor;

  // // Compute deviatoric trial stress,  sigma_tr_dev
  // Matrix<Real> sigma_tr_dev(sigma_tr);
  // sigma_tr_dev -= Matrix<Real>::eye(dim, sigma_tr.trace() / 3.0);

  // // Compute effective deviatoric trial stress
  // Real s = sigma_tr_dev.doubleDot(sigma_tr_dev);
  // Real sigma_tr_dev_eff=std::sqrt(3./2. * s);

  // // Compute deviatoric stress,  sigma_dev
  // Matrix<Real> sigma_dev(sigma_tensor);
  // sigma_dev -= Matrix<Real>::eye(dim, sigma_tensor.trace() / 3.0);

  // // Compute effective deviatoric stress
  // s =  sigma_dev.doubleDot(sigma_dev);
  // Real sigma_dev_eff = std::sqrt(3./2. * s);

  // Real xr = 0.0;
  // if(sigma_tr_dev_eff > sigma_dev_eff * Math::getTolerance())
  //   xr = sigma_dev_eff / sigma_tr_dev_eff;

  // Real __attribute__((unused)) q = 1.5 * (1. / (1. +  3. * this->mu  / this->h) - xr);

  UInt cols = tangent.cols();
  UInt rows = tangent.rows();

  for (UInt m = 0; m < rows; ++m) {
    UInt i = VoigtHelper<dim>::vec[m][0];
    UInt j = VoigtHelper<dim>::vec[m][1];

    for (UInt n = 0; n < cols; ++n) {
      UInt k = VoigtHelper<dim>::vec[n][0];
      UInt l = VoigtHelper<dim>::vec[n][1];

      // This section of the code is commented
      // There were some problems with the convergence of plastic-coupled simulations with thermal expansion
      // XXX: DO NOT REMOVE
      /*if (((sigma_tr_dev_eff-iso_hardening-sigmay) > 0) && (xr > 0)) {
        tangent(m,n) =
        2. * this->mu * q * (sigma_tr_dev (i,j) / sigma_tr_dev_eff) * (sigma_tr_dev (k,l) / sigma_tr_dev_eff) +
        (i==k) * (j==l) * 2. * this->mu * xr +
        (i==j) * (k==l) * (this->kpa - 2./3. * this->mu * xr);
        if ((m == n) && (m>=dim))
        tangent(m, n) = tangent(m, n) - this->mu * xr;
        } else {*/
      tangent(m,n) =  (i==k) * (j==l) * 2. * this->mu +
        (i==j) * (k==l) * this->lambda;
      tangent(m,n) -= (m==n) * (m>=dim) * this->mu;
      //}
      //correct tangent stiffness for shear component
    }
  }
}

