/**
 * @file   remove_damaged_with_damage_rate_weight_function.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Cyprien Wolff <cyprien.wolff@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  Removed damaged weight function for non local materials
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
#include "base_weight_function.hh"
/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_REMOVE_DAMAGED_WITH_DAMAGE_RATE_WEIGHT_FUNCTION_HH__
#define __AKANTU_REMOVE_DAMAGED_WITH_DAMAGE_RATE_WEIGHT_FUNCTION_HH__

namespace akantu {
/* -------------------------------------------------------------------------- */
/* Remove damaged with damage rate weight function */
/* -------------------------------------------------------------------------- */

class RemoveDamagedWithDamageRateWeightFunction : public BaseWeightFunction {
public:
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
  RemoveDamagedWithDamageRateWeightFunction(NonLocalManager & manager)
      : BaseWeightFunction(manager, "remove_damage_with_damage_rate"),
        damage_with_damage_rate(nullptr) {
    this->registerParam<Real>("damage_limit",
                              this->damage_limit_with_damage_rate, 1,
                              _pat_parsable, "Damage Threshold");
    this->init();
  }

  /* --------------------------------------------------------------------------
   */
  /* Base Weight Function inherited methods */
  /* --------------------------------------------------------------------------
   */
  inline Real operator()(Real r,
                         const __attribute__((unused)) IntegrationPoint & q1,
                         const IntegrationPoint & q2);

  inline void init() override;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// limit at which a point is considered as complitely broken
  Real damage_limit_with_damage_rate;

  /// internal pointer to the current damage vector
  ElementTypeMapReal * damage_with_damage_rate;
};

} // namespace akantu

#include "remove_damaged_with_damage_rate_weight_function_inline_impl.hh"

#endif /* __AKANTU_REMOVE_DAMAGED_WITH_DAMAGE_WEIGHT_FUNCTION_HH__ */
