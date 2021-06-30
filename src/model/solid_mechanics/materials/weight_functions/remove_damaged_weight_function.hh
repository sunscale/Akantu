/**
 * @file   remove_damaged_weight_function.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Cyprien Wolff <cyprien.wolff@epfl.ch>
 *
 * @date creation: Mon Aug 24 2015
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  Removed damaged weight function for non local materials
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
#include "base_weight_function.hh"
/* -------------------------------------------------------------------------- */
#ifndef AKANTU_REMOVE_DAMAGED_WEIGHT_FUNCTION_HH_
#define AKANTU_REMOVE_DAMAGED_WEIGHT_FUNCTION_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
/*  Remove damaged weight function                                            */
/* -------------------------------------------------------------------------- */
class RemoveDamagedWeightFunction : public BaseWeightFunction {
public:
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
  RemoveDamagedWeightFunction(NonLocalManager & manager)
      : BaseWeightFunction(manager, "remove_damaged"), damage(nullptr) {
    this->registerParam("damage_limit", this->damage_limit, 1., _pat_parsable,
                        "Damage Threshold");
    this->init();
  }

  /* --------------------------------------------------------------------------
   */
  /* Base Weight Function inherited methods */
  /* --------------------------------------------------------------------------
   */

  inline void init() override;

  inline Real operator()(Real r, const IntegrationPoint & q1,
                         const IntegrationPoint & q2);

  /* ------------------------------------------------------------------------ */
  /* Data Accessor inherited members                                          */
  /* ------------------------------------------------------------------------ */

  inline UInt getNbData(const Array<Element> & elements,
                        const SynchronizationTag & tag) const override;

  inline void packData(CommunicationBuffer & buffer,
                       const Array<Element> & elements,
                       const SynchronizationTag & tag) const override;

  inline void unpackData(CommunicationBuffer & buffer,
                         const Array<Element> & elements,
                         const SynchronizationTag & tag) override;

private:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
  /// limit at which a point is considered as complitely broken
  Real damage_limit;

  /// internal pointer to the current damage vector
  ElementTypeMapReal * damage;
};

} // namespace akantu

#include "remove_damaged_weight_function_inline_impl.hh"

#endif /* AKANTU_REMOVE_DAMAGED_WEIGHT_FUNCTION_HH_ */
