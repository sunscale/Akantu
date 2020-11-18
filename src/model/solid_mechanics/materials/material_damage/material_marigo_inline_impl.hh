/**
 * @file   material_marigo_inline_impl.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Aug 04 2010
 * @date last modification: Fri Dec 02 2016
 *
 * @brief  Implementation of the inline functions of the material marigo
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_marigo.hh"

#ifndef AKANTU_MATERIAL_MARIGO_INLINE_IMPL_HH_
#define AKANTU_MATERIAL_MARIGO_INLINE_IMPL_HH_

namespace akantu {

template <UInt spatial_dimension>
inline void MaterialMarigo<spatial_dimension>::computeStressOnQuad(
    Matrix<Real> & grad_u, Matrix<Real> & sigma, Real & dam, Real & Y,
    Real & Ydq) {
  MaterialElastic<spatial_dimension>::computeStressOnQuad(grad_u, sigma);

  Y = 0;
  for (UInt i = 0; i < spatial_dimension; ++i) {
    for (UInt j = 0; j < spatial_dimension; ++j) {
      Y += sigma(i, j) * (grad_u(i, j) + grad_u(j, i)) / 2.;
    }
  }
  Y *= 0.5;

  if (damage_in_y) {
    Y *= (1 - dam);
  }

  if (yc_limit) {
    Y = std::min(Y, Yc);
  }

  if (!this->is_non_local) {
    computeDamageAndStressOnQuad(sigma, dam, Y, Ydq);
  }
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
inline void MaterialMarigo<spatial_dimension>::computeDamageAndStressOnQuad(
    Matrix<Real> & sigma, Real & dam, Real & Y, Real & Ydq) {
  Real Fd = Y - Ydq - Sd * dam;

  if (Fd > 0) {
    dam = (Y - Ydq) / Sd;
  }
  dam = std::min(dam, Real(1.));

  sigma *= 1 - dam;
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
inline UInt MaterialMarigo<spatial_dimension>::getNbData(
    const Array<Element> & elements, const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  if (tag == SynchronizationTag::_smm_init_mat) {
    size += sizeof(Real) * this->getModel().getNbIntegrationPoints(elements);
  }

  size += MaterialDamage<spatial_dimension>::getNbData(elements, tag);

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
inline void MaterialMarigo<spatial_dimension>::packData(
    CommunicationBuffer & buffer, const Array<Element> & elements,
    const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  if (tag == SynchronizationTag::_smm_init_mat) {
    this->packElementDataHelper(Yd, buffer, elements);
  }

  MaterialDamage<spatial_dimension>::packData(buffer, elements, tag);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
inline void
MaterialMarigo<spatial_dimension>::unpackData(CommunicationBuffer & buffer,
                                              const Array<Element> & elements,
                                              const SynchronizationTag & tag) {
  AKANTU_DEBUG_IN();

  if (tag == SynchronizationTag::_smm_init_mat) {
    this->unpackElementDataHelper(Yd, buffer, elements);
  }

  MaterialDamage<spatial_dimension>::unpackData(buffer, elements, tag);

  AKANTU_DEBUG_OUT();
}

} // namespace akantu

#endif /* AKANTU_MATERIAL_MARIGO_INLINE_IMPL_HH_ */
