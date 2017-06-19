/**
 * @file   damaged_weight_function.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Cyprien Wolff <cyprien.wolff@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Thu Oct 15 2015
 *
 * @brief  Damaged weight function for non local materials
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "base_weight_function.hh"
/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_DAMAGED_WEIGHT_FUNCTION_HH__
#define __AKANTU_DAMAGED_WEIGHT_FUNCTION_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
/*  Damage weight function                                                    */
/* -------------------------------------------------------------------------- */
class DamagedWeightFunction : public BaseWeightFunction {
public:
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
  DamagedWeightFunction(NonLocalManager & manager) : BaseWeightFunction(manager, "damaged"),
						     damage(NULL) {
    this->init();
  }

  /* -------------------------------------------------------------------------- */
  /* Base Weight Function inherited methods                                     */
  /* -------------------------------------------------------------------------- */

  /// set the pointers of internals to the right flattend version
  virtual void init();

  inline Real operator()(Real r,
			 const __attribute__((unused)) IntegrationPoint & q1,
			 const IntegrationPoint & q2);

private:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

  /// internal pointer to the current damage vector
  ElementTypeMapReal * damage;

};


#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "damaged_weight_function_inline_impl.cc"
#endif

} // akantu

#endif /* __AKANTU_DAMAGED_WEIGHT_FUNCTION_HH__ */
