/**
 * @file   damaged_weight_function_inline_impl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Cyprien Wolff <cyprien.wolff@epfl.ch>
 *
 * @date creation: Mon Aug 24 2015
 * @date last modification: Wed Feb 03 2016
 *
 * @brief  Implementation of inline function of damaged weight function
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
#include "damaged_weight_function.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
inline Real DamagedWeightFunction::operator()(Real r,
                                              const __attribute__((unused))
                                              IntegrationPoint & q1,
                                              const IntegrationPoint & q2) {
  /// compute the weight
  UInt quad = q2.global_num;
  Array<Real> & dam_array = (*this->damage)(q2.type, q2.ghost_type);
  Real D = dam_array(quad);
  Real Radius_t = 0;
  Real Radius_init = this->R2;

  //    if(D <= 0.5)
  //      {
  //	Radius_t = 2*D*Radius_init;
  //      }
  //    else
  //      {
  //	Radius_t = 2*Radius_init*(1-D);
  //      }
  //

  Radius_t = Radius_init * (1 - D);

  Radius_init *= Radius_init;
  Radius_t *= Radius_t;

  if (Radius_t < Math::getTolerance()) {
    Radius_t = 0.001 * Radius_init;
  }

  Real expb =
      (2 * std::log(0.51)) / (std::log(1.0 - 0.49 * Radius_t / Radius_init));
  Int expb_floor = std::floor(expb);
  Real b = expb_floor + expb_floor % 2;
  Real alpha = std::max(0., 1. - r * r / Radius_init);
  Real w = std::pow(alpha, b);
  return w;
}

/* -------------------------------------------------------------------------- */
inline void DamagedWeightFunction::init() {
  this->damage = &(this->manager.registerWeightFunctionInternal("damage"));
}
} // namespace akantu
