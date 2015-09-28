/**
 * @file   non_local_neighborhood.hh
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Sat Sep 26 18:29:59 2015
 *
 * @brief Non-local neighborhood for non-local averaging based on
 * weight function
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
#ifndef __AKANTU_NON_LOCAL_NEIGHBORHOOD_HH__
#define __AKANTU_NON_LOCAL_NEIGHBORHOOD_HH__
/* -------------------------------------------------------------------------- */
#include "base_weight_function.hh"
#include "non_local_neighborhood_base.hh"
/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__

class NonLocalNeighborhood : NonLocalNeighborhoodBase {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  NonLocalNeighborhood(SolidMechanicsModel & model, Real radius, const ID & weight_type);
  virtual ~NonLocalNeighborhood();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  
  void computeWeights();

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// the weights associated to the pairs
  Array<Real> * pair_weight[2];

  /// count the number of calls of computeStress
  UInt compute_stress_calls;

  /// weight function
  BaseWeightFunction * weight_function;

};


__END_AKANTU__

#endif /* __AKANTU_NON_LOCAL_NEIGHBORHOOD_HH__ */
