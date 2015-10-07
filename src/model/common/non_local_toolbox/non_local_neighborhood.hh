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
#include "parsable.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
  class NonLocalManager;
  class TestWeightFunction;
}

__BEGIN_AKANTU__


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


class TestWeightFunction : public Parsable {

public:
  TestWeightFunction(NonLocalManager & manager) : Parsable(_st_weight_function, "weight_function:"), 
							manager(manager)
{
    this->registerParam("update_rate"  , update_rate, 0U  ,
			_pat_parsmod, "Update frequency");

  }

  void setRadius(Real radius) {
    /// set the non-local radius and update R^2 accordingly
    this->R = radius; 
    this->R2 = this->R * this->R;
  }

Real operator()(Real r,
		const __attribute__((unused)) QuadraturePoint & q1,
		const __attribute__((unused)) QuadraturePoint & q2) {
  /// initialize the weightx
  Real w = 0;
  /// compute weight for given r 
  if(r <= this->R) {
    Real alpha = (1. - r*r / this->R2);
    w = alpha * alpha;
    // *weight = 1 - sqrt(r / radius);
  }
  return w;
}
  void updateInternals(){};
  UInt getUpdateRate(){return update_rate;};

protected:

  NonLocalManager & manager;
  Real R;
  Real R2;
  UInt update_rate;

};

// class TestDamageWeightFunction : public TestWeightFunction {

// public:
//   TestDamageWeightFunction(NonLocalManager & manager) : TestWeightFunction(manager){ 
//     this->registerParam("damage_limit", this->damage_limit, 1., _pat_parsable, "Damage Threshold");
//     this->setInternals();}

// Real operator()(Real r,
// 		const __attribute__((unused)) QuadraturePoint & q1,
// 		const __attribute__((unused)) QuadraturePoint & q2) {
//   /// compute the weight
//   UInt quad = q2.global_num;

//   if(q1 == q2) return 1.;

//   Array<Real> & dam = (*(this->damage))(q2.type, q2.ghost_type);
//   Real D = dam(quad);
//   Real w = 0.;
//   if(D < damage_limit) {
//     Real alpha = std::max(0., 1. - r*r / this->R2);
//     w = alpha * alpha;
//   }
//   return w;
// }
//   void updateInternals(){};
//   UInt getUpdateRate(){return update_rate;};

//   void setInternals() {
//     this->damage = &(this->manager.registerWeightFunctionInternal("damage"));}

// protected:

//  /// limit at which a point is considered as complitely broken
//   Real damage_limit;

//   /// internal pointer to the current damage vector
//   ElementTypeMapReal * damage;

// };


__END_AKANTU__

/* -------------------------------------------------------------------------- */
/* Implementation of template functions                                       */
/* -------------------------------------------------------------------------- */
#include "non_local_neighborhood_tmpl.hh"
/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "non_local_neighborhood_inline_impl.cc"



#endif /* __AKANTU_NON_LOCAL_NEIGHBORHOOD_HH__ */

