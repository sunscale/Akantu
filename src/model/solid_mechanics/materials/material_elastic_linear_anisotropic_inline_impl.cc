/**
 * @file   material_elastic_linear_anisotropic_inline_impl.cc
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
 * @section LICENSE
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

#ifndef __AKANTU_MATERIAL_ELASTIC_LINEAR_ANISOTROPIC_INLINE_IMPL_CC__
#define __AKANTU_MATERIAL_ELASTIC_LINEAR_ANISOTROPIC_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt dim>
inline void MaterialElasticLinearAnisotropic<dim>::computeStressOnQuad(
    const Matrix<Real> & grad_u, Matrix<Real> & sigma) const {
  // Wikipedia convention:
  // 2*eps_ij (i!=j) = voigt_eps_I
  // http://en.wikipedia.org/wiki/Voigt_notation
  Vector<Real> voigt_strain(voigt_h::size);
  Vector<Real> voigt_stress(voigt_h::size);

  for (UInt I = 0; I < voigt_h::size; ++I) {
    Real voigt_factor = voigt_h::factors[I];
    UInt i = voigt_h::vec[I][0];
    UInt j = voigt_h::vec[I][1];

    voigt_strain(I) = voigt_factor * (grad_u(i, j) + grad_u(j, i)) / 2.;
  }

  voigt_stress = this->C * voigt_strain;

  for (UInt I = 0; I < voigt_h::size; ++I) {
    UInt i = voigt_h::vec[I][0];
    UInt j = voigt_h::vec[I][1];

    sigma(i, j) = sigma(j, i) = voigt_stress(I);
  }
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

#endif /* __AKANTU_MATERIAL_ELASTIC_LINEAR_ANISOTROPIC_INLINE_IMPL_CC__ */
