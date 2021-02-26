/**
 * @file   base_weight_function_inline_impl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Cyprien Wolff <cyprien.wolff@epfl.ch>
 *
 * @date creation: Wed Sep 01 2010
 * @date last modification: Wed Sep 27 2017
 *
 * @brief  Implementation of inline function of base weight function
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
#include "base_weight_function.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_BASE_WEIGHT_FUNCTION_INLINE_IMPL_HH_
#define AKANTU_BASE_WEIGHT_FUNCTION_INLINE_IMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
inline void BaseWeightFunction::init() {
  /// compute R^2 for a given non-local radius
  this->R2 = this->R * this->R;
}

/* -------------------------------------------------------------------------- */
inline void BaseWeightFunction::setRadius(Real radius) {
  /// set the non-local radius and update R^2 accordingly
  this->R = radius;
  this->R2 = this->R * this->R;
}

/* -------------------------------------------------------------------------- */
inline Real
BaseWeightFunction::operator()(Real r, const IntegrationPoint & /* q1 */,
                               const IntegrationPoint & /* q2 */) const {

  /// initialize the weight
  Real w = 0;
  /// compute weight for given r
  if (r <= this->R) {
    Real alpha = (1. - r * r / this->R2);
    w = alpha * alpha;
    // *weight = 1 - sqrt(r / radius);
  }

  return w;
}

} // namespace akantu
#endif /* AKANTU_BASE_WEIGHT_FUNCTION_INLINE_IMPL_HH_ */
