/**
 * @file   material_plastic_inline_impl.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Apr 07 2014
 * @date last modification: Wed Jan 24 2018
 *
 * @brief  Implementation of the inline functions of akantu::MaterialPlastic
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#ifndef MATERIAL_PLASTIC_INLINE_IMPL_H
#define MATERIAL_PLASTIC_INLINE_IMPL_H

#include "material_plastic.hh"

namespace akantu {

template <UInt dim>
inline void MaterialPlastic<dim>::computeStressAndInelasticStrainOnQuad(
    const Matrix<Real> & delta_grad_u, Matrix<Real> & sigma,
    const Matrix<Real> & previous_sigma, Matrix<Real> & inelastic_strain,
    const Matrix<Real> & previous_inelastic_strain,
    const Matrix<Real> & delta_inelastic_strain) const {
  Matrix<Real> grad_u_elastic(dim, dim);
  grad_u_elastic.copy(delta_grad_u);
  grad_u_elastic -= delta_inelastic_strain;
  Matrix<Real> sigma_elastic(dim, dim);
  MaterialElastic<dim>::computeStressOnQuad(grad_u_elastic, sigma_elastic);
  sigma.copy(previous_sigma);
  sigma += sigma_elastic;
  inelastic_strain.copy(previous_inelastic_strain);
  inelastic_strain += delta_inelastic_strain;
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
inline void MaterialPlastic<dim>::computeStressAndInelasticStrainOnQuad(
    const Matrix<Real> & grad_u, const Matrix<Real> & previous_grad_u,
    Matrix<Real> & sigma, const Matrix<Real> & previous_sigma,
    Matrix<Real> & inelastic_strain,
    const Matrix<Real> & previous_inelastic_strain,
    const Matrix<Real> & delta_inelastic_strain) const {
  Matrix<Real> delta_grad_u(grad_u);
  delta_grad_u -= previous_grad_u;

  computeStressAndInelasticStrainOnQuad(
      delta_grad_u, sigma, previous_sigma, inelastic_strain,
      previous_inelastic_strain, delta_inelastic_strain);
}

} // namespace akantu

#endif /* MATERIAL_PLASTIC_INLINE_IMPL_H */
