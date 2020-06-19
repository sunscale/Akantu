/**
 * @file   remove_damaged_weight_function_inline_impl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Cyprien Wolff <cyprien.wolff@epfl.ch>
 *
 * @date creation: Mon Aug 24 2015
 * @date last modification: Thu Jul 06 2017
 *
 * @brief  Implementation of inline function of remove damaged weight function
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "remove_damaged_weight_function.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_REMOVE_DAMAGED_WEIGHT_FUNCTION_INLINE_IMPL_HH__
#define __AKANTU_REMOVE_DAMAGED_WEIGHT_FUNCTION_INLINE_IMPL_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
inline Real RemoveDamagedWeightFunction::
operator()(Real r, const __attribute__((unused)) IntegrationPoint & q1,
           const IntegrationPoint & q2) {
  /// compute the weight
  UInt quad = q2.global_num;

  if (q1 == q2)
    return 1.;

  Array<Real> & dam_array = (*this->damage)(q2.type, q2.ghost_type);
  Real D = dam_array(quad);
  Real w = 0.;
  if (D < damage_limit * (1 - Math::getTolerance())) {
    Real alpha = std::max(0., 1. - r * r / this->R2);
    w = alpha * alpha;
  }
  return w;
}

/* -------------------------------------------------------------------------- */
inline void RemoveDamagedWeightFunction::init() {
  this->damage = &(this->manager.registerWeightFunctionInternal("damage"));
}

/* -------------------------------------------------------------------------- */
inline UInt
RemoveDamagedWeightFunction::getNbData(const Array<Element> & elements,
                                       const SynchronizationTag & tag) const {

  if (tag == SynchronizationTag::_mnl_weight)
    return this->manager.getModel().getNbIntegrationPoints(elements) *
           sizeof(Real);

  return 0;
}

/* -------------------------------------------------------------------------- */
inline void
RemoveDamagedWeightFunction::packData(CommunicationBuffer & buffer,
                                      const Array<Element> & elements,
                                      const SynchronizationTag & tag) const {
  if (tag == SynchronizationTag::_mnl_weight) {
    DataAccessor<Element>::packElementalDataHelper<Real>(
        *damage, buffer, elements, true,
        this->manager.getModel().getFEEngine());
  }
}

/* -------------------------------------------------------------------------- */
inline void
RemoveDamagedWeightFunction::unpackData(CommunicationBuffer & buffer,
                                        const Array<Element> & elements,
                                        const SynchronizationTag & tag) {
  if (tag == SynchronizationTag::_mnl_weight) {
    DataAccessor<Element>::unpackElementalDataHelper<Real>(
        *damage, buffer, elements, true,
        this->manager.getModel().getFEEngine());
  }
}

} // namespace akantu

#endif /* __AKANTU_REMOVE_DAMAGED_WEIGHT_FUNCTION_INLINE_IMPL_HH__ */
