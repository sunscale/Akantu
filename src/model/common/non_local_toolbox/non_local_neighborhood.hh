/**
 * @file   non_local_neighborhood.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Nov 25 2015
 *
 * @brief  Non-local neighborhood for non-local averaging based on
 * weight function
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
#ifndef __AKANTU_NON_LOCAL_NEIGHBORHOOD_HH__
#define __AKANTU_NON_LOCAL_NEIGHBORHOOD_HH__
/* -------------------------------------------------------------------------- */
#include "base_weight_function.hh"
#include "non_local_neighborhood_base.hh"
#include "parsable.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
  class NonLocalManager;
  class TestWeightFunction;
}

namespace akantu {


template<class WeightFunction = BaseWeightFunction>
class NonLocalNeighborhood : public NonLocalNeighborhoodBase {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  NonLocalNeighborhood(NonLocalManager & manager, 
		       const ElementTypeMapReal & quad_coordinates,
		       const ID & id = "neighborhood",
		       const MemoryID & memory_id = 0);
  virtual ~NonLocalNeighborhood();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

   /// compute the weights for non-local averaging
  void computeWeights();
  
  /// save the pair of weights in a file
  void saveWeights(const std::string & filename) const;

  /// compute the non-local counter part for a given element type map
  virtual void weightedAverageOnNeighbours(const ElementTypeMapReal & to_accumulate,
					   ElementTypeMapReal & accumulated,
					   UInt nb_degree_of_freedom,
					   const GhostType & ghost_type2) const;

  /// update the weights based on the weight function
  void updateWeights();

  
  /// register a new non-local variable in the neighborhood
  virtual void registerNonLocalVariable(const ID & id);

protected:
  virtual inline UInt getNbDataForElements(const Array<Element> & elements,
					  SynchronizationTag tag) const;

  virtual inline void packElementData(CommunicationBuffer & buffer,
				      const Array<Element> & elements,
				      SynchronizationTag tag) const;

  virtual inline void unpackElementData(CommunicationBuffer & buffer,
					const Array<Element> & elements,
					SynchronizationTag tag);

  /* -------------------------------------------------------------------------- */
  /* Accessor                                                                   */
  /* -------------------------------------------------------------------------- */
  AKANTU_GET_MACRO(NonLocalManager, *non_local_manager, const NonLocalManager &);
  AKANTU_GET_MACRO_NOT_CONST(NonLocalManager, *non_local_manager, NonLocalManager &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// Pointer to non-local manager class
  NonLocalManager * non_local_manager;

  /// the weights associated to the pairs
  Array<Real> * pair_weight[2];

  /// weight function
  WeightFunction * weight_function;

};

} // akantu

/* -------------------------------------------------------------------------- */
/* Implementation of template functions                                       */
/* -------------------------------------------------------------------------- */
#include "non_local_neighborhood_tmpl.hh"
/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "non_local_neighborhood_inline_impl.cc"



#endif /* __AKANTU_NON_LOCAL_NEIGHBORHOOD_HH__ */

