/**
 * @file   remove_damaged_weight_function_inline_impl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Cyprien Wolff <cyprien.wolff@epfl.ch>
 *
 * @date creation: Fri Apr 13 2012
 * @date last modification: Thu Jun 05 2014
 *
 * @brief Implementation of inline function of remove damaged weight function 
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
inline void RemoveDamagedWeightFunction::selectType(__attribute__((unused)) ElementType type1,
						    __attribute__((unused)) GhostType ghost_type1,
						    ElementType type2,
						    GhostType ghost_type2) {
  /// select the damage array for a given type: For optimization
  selected_damage = &(this->material.template getArray<Real>("damage", type2, ghost_type2));
  //    selected_damage = &mat.getDamage(type2, ghost_type2);
}

/* -------------------------------------------------------------------------- */
inline Real RemoveDamagedWeightFunction::operator()(Real r,
						    const __attribute__((unused)) QuadraturePoint & q1,
						    const QuadraturePoint & q2) {
  /// compute the weight
  UInt quad = q2.global_num;

  if(q1 == q2) return 1.;

  Real D = (*selected_damage)(quad);
  Real w = 0.;
  if(D < damage_limit) {
    Real alpha = std::max(0., 1. - r*r / this->R2);
    w = alpha * alpha;
  }
  return w;
}

/* -------------------------------------------------------------------------- */
inline UInt RemoveDamagedWeightFunction::getNbDataForElements(const Array<Element> & elements,
							      SynchronizationTag tag) const {
  if(tag == _gst_mnl_weight)
    return this->material.getModel().getNbQuadraturePoints(elements) * sizeof(Real);
  
  return 0;
}

/* -------------------------------------------------------------------------- */
inline void RemoveDamagedWeightFunction::packElementData(CommunicationBuffer & buffer,
							 const Array<Element> & elements,
							 SynchronizationTag tag) const {
  if(tag == _gst_mnl_weight) {
    ElementTypeMapArray<Real> & damage = this->material.template getInternal<Real>("damage");
    this->material.packElementDataHelper(damage,
					 buffer,
					 elements);
  }
}

/* -------------------------------------------------------------------------- */
inline void RemoveDamagedWeightFunction::unpackElementData(CommunicationBuffer & buffer,
							   const Array<Element> & elements,
							   SynchronizationTag tag) {
  if(tag == _gst_mnl_weight) {
    ElementTypeMapArray<Real> & damage = this->material.template getInternal<Real>("damage");
    this->material.unpackElementDataHelper(damage,
					   buffer,
					   elements);
  }
}




