/**
 * @file   local_material_damage_inline_impl.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Aug 10 2015
 * @date last modification: Mon Jan 18 2016
 *
 * @brief  Implementation of the inline functions of the material damage
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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

/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_LOCAL_MATERIAL_DAMAGE_INLINE_IMPL_HH__
#define __AKANTU_LOCAL_MATERIAL_DAMAGE_INLINE_IMPL_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
inline void LocalMaterialDamage::computeStressOnQuad(Matrix<Real> & grad_u,
                                                     Matrix<Real> & sigma,
                                                     Real & dam) {

  Real trace = grad_u.trace();

  /// \sigma_{ij} = \lambda * (\nabla u)_{kk} * \delta_{ij} + \mu * (\nabla
  /// u_{ij} + \nabla u_{ji})
  for (UInt i = 0; i < spatial_dimension; ++i) {
    for (UInt j = 0; j < spatial_dimension; ++j) {
      sigma(i, j) =
          (i == j) * lambda * trace + mu * (grad_u(i, j) + grad_u(j, i));
    }
  }

  Real Y = 0;
  for (UInt i = 0; i < spatial_dimension; ++i) {
    for (UInt j = 0; j < spatial_dimension; ++j) {
      Y += sigma(i, j) * grad_u(i, j);
    }
  }
  Y *= 0.5;

  Real Fd = Y - Yd - Sd * dam;

  if (Fd > 0)
    dam = (Y - Yd) / Sd;
  dam = std::min(dam, 1.);

  sigma *= 1 - dam;
}

/* -------------------------------------------------------------------------- */
inline void LocalMaterialDamage::computePotentialEnergyOnQuad(
    Matrix<Real> & grad_u, Matrix<Real> & sigma, Real & epot) {
  epot = 0.;
  for (UInt i = 0, t = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j, ++t)
      epot += sigma(i, j) * (grad_u(i, j) - (i == j));
  epot *= .5;
}

/* -------------------------------------------------------------------------- */
inline Real LocalMaterialDamage::getCelerity(__attribute__((unused))
                                             const Element & element) const {
  // Here the fastest celerity is the push wave speed
  return (std::sqrt((2 * mu + lambda) / rho));
}

} // namespace akantu

#endif /* __AKANTU_LOCAL_MATERIAL_DAMAGE_INLINE_IMPL_HH__ */
