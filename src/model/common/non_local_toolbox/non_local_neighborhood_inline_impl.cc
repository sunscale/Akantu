/**
 * @file   non_local_neighborhood_inline_impl.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Mon Oct  5 09:36:33 2015
 *
 * @brief  Implementation of inline functions of non-local neighborhood class
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
#ifndef __AKANTU_NON_LOCAL_NEIGHBORHOOD_INLINE_IMPL_CC__
#define __AKANTU_NON_LOCAL_NEIGHBORHOOD_INLINE_IMPL_CC__

__BEGIN_AKANTU__
/* -------------------------------------------------------------------------- */
template<class WeightFunction>
inline UInt NonLocalNeighborhood<WeightFunction>::getNbDataForElements(const Array<Element> & elements,
								       SynchronizationTag tag) const {
  UInt size = 0;

  if(tag == _gst_mnl_for_average) {
    std::set<ID>::const_iterator it = non_local_variables.begin();
    std::set<ID>::const_iterator end = non_local_variables.end();
   
    for(;it != end; ++it) {   
      size += this->non_local_manager->getNbDataForElements(elements, *it);
    }
  }

  size += this->weight_function->getNbDataForElements(elements, tag);

  return size;
}

/* -------------------------------------------------------------------------- */
template<class WeightFunction>
inline void NonLocalNeighborhood<WeightFunction>::packElementData(CommunicationBuffer & buffer,
								  const Array<Element> & elements,
								  SynchronizationTag tag) const {
  if(tag == _gst_mnl_for_average) {

    std::set<ID>::const_iterator it = non_local_variables.begin();
    std::set<ID>::const_iterator end = non_local_variables.end();

    for(;it != end; ++it) {
      this->non_local_manager->packElementData(buffer, elements, tag, *it);
    }
  }

  this->weight_function->packElementData(buffer, elements, tag);
}

/* -------------------------------------------------------------------------- */
template<class WeightFunction>
inline void NonLocalNeighborhood<WeightFunction>::unpackElementData(CommunicationBuffer & buffer,
								    const Array<Element> & elements,
								    SynchronizationTag tag) {
  if(tag == _gst_mnl_for_average) {

    std::set<ID>::const_iterator it = non_local_variables.begin();
    std::set<ID>::const_iterator end = non_local_variables.end();

    for(;it != end; ++it) {
      this->non_local_manager->unpackElementData(buffer, elements, tag, *it);
    }
  }
 
  this->weight_function->unpackElementData(buffer, elements, tag);
}

__END_AKANTU__

#endif /* __AKANTU_NON_LOCAL_NEIGHBORHOOD_INLINE_IMPL_CC__ */
