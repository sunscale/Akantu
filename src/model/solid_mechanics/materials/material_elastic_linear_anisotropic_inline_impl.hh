/**
 * @file   material_elastic_linear_anisotropic_inline_impl.hh
 *
 * @author Enrico Milanese <enrico.milanese@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Feb 16 2018
 * @date last modification: Fri Feb 16 2018
 *
 * @brief  Implementation of the inline functions of the material elastic linear
 * anisotropic
 *
 *
 * Copyright (©) 2016-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_elastic_linear_anisotropic.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_ELASTIC_LINEAR_ANISOTROPIC_INLINE_IMPL_HH_
#define AKANTU_MATERIAL_ELASTIC_LINEAR_ANISOTROPIC_INLINE_IMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt dim>
inline void MaterialElasticLinearAnisotropic<dim>::computeStressOnQuad(
    const Matrix<Real> & grad_u, Matrix<Real> & sigma) const {
  auto voigt_strain = strainToVoigt<dim>(gradUToEpsilon<dim>(grad_u));
  auto voigt_stress = this->C * voigt_strain;
  voigtToStress<dim>(voigt_stress, sigma);
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
inline void MaterialElasticLinearAnisotropic<dim>::computeTangentModuliOnQuad(
    Matrix<Real> & tangent) const {

  tangent.copy(this->C);
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
inline void MaterialElasticLinearAnisotropic<dim>::computePotentialEnergyOnQuad(
    const Matrix<Real> & grad_u, const Matrix<Real> & sigma, Real & epot) {

  AKANTU_DEBUG_ASSERT(this->symmetric,
                      "The elastic constants matrix is not symmetric,"
                      "energy is not path independent.");

  epot = .5 * sigma.doubleDot(grad_u);
}

} // namespace akantu

#endif /* AKANTU_MATERIAL_ELASTIC_LINEAR_ANISOTROPIC_INLINE_IMPL_HH_ */
