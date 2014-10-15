/**
 * @file   material_viscoplastic.cc
 *
 * @author Ramin Aghababaei <ramin.aghababaei@epfl.ch>
 *
 * @date   Tue Jul 09 18:15:37 20130
 *
 * @brief  Implementation of the inline functions of the material viscoplastic
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#include <cmath>
#include "material_viscoplastic.hh"


/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template<UInt dim>
inline void MaterialViscoPlastic<dim>::computeStressOnQuad(const Matrix<Real> & grad_u,
                                                           const Matrix<Real> & previous_grad_u,
                                                           Matrix<Real> & sigma,
                                                           const Matrix<Real> & previous_sigma,
                                                           Matrix<Real> & inelastic_strain,
                                                           const Matrix<Real> & previous_inelastic_strain,
                                                           Real & iso_hardening) const {
  // Infinitesimal plastic
  // Compute stress magnitude
  Real s = sigma.doubleDot(sigma);
  Real sigma_mag=sqrt(s);

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

  //Update resistance stress
  iso_hardening = iso_hardening + this->h * dep_mag;

  MaterialPlastic<dim>::computeStressAndInelasticStrainOnQuad(grad_delta_u, sigma, previous_sigma,
                                                              inelastic_strain, previous_inelastic_strain,
                                                              delta_inelastic_strain);

}

/* -------------------------------------------------------------------------- */
template<UInt dim>
inline void MaterialViscoPlastic<dim>::computeTangentModuliOnQuad(Matrix<Real> & tangent,
                                                                  const Matrix<Real> & grad_u,
                                                                  const Matrix<Real> & previous_grad_u,
                                                                  const Matrix<Real> & sigma_tensor,
                                                                  const Matrix<Real> & previous_sigma_tensor,
                                                                  const Real & iso_hardening) const {
  UInt cols = tangent.cols();
  UInt rows = tangent.rows();

  for (UInt m = 0; m < rows; ++m) {
    UInt i = VoigtHelper<dim>::vec[m][0];
    UInt j = VoigtHelper<dim>::vec[m][1];

    for (UInt n = 0; n < cols; ++n) {
      UInt k = VoigtHelper<dim>::vec[n][0];
      UInt l = VoigtHelper<dim>::vec[n][1];
      tangent(m,n) = (i==k) * (j==l) * 2. * this->mu + (i==j) * (k==l) * this->lambda;
      tangent(m,n) -= (m==n) * (m>=dim) * this->mu;
    }
  }
}
