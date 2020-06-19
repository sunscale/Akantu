/**
 * @file   neighborhood_max_criterion_inline_impl.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sat Sep 26 2015
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  Implementation of inline functions for class NeighborhoodMaxCriterion
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
#include "model.hh"
#include "neighborhood_max_criterion.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_NEIGHBORHOOD_MAX_CRITERION_INLINE_IMPL_HH__
#define __AKANTU_NEIGHBORHOOD_MAX_CRITERION_INLINE_IMPL_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
inline UInt
NeighborhoodMaxCriterion::getNbDataForElements(const Array<Element> & elements,
                                               SynchronizationTag tag) const {
  UInt nb_quadrature_points = this->model.getNbIntegrationPoints(elements);
  UInt size = 0;

  if (tag == SynchronizationTag::_nh_criterion) {
    size += sizeof(Real) * nb_quadrature_points;
  }

  return size;
}

/* -------------------------------------------------------------------------- */
inline void
NeighborhoodMaxCriterion::packElementData(CommunicationBuffer & buffer,
                                          const Array<Element> & elements,
                                          SynchronizationTag tag) const {
  if (tag == SynchronizationTag::_nh_criterion) {
    this->packElementalDataHelper(criterion, buffer, elements, true,
                                  this->model.getFEEngine());
  }
}

/* -------------------------------------------------------------------------- */
inline void
NeighborhoodMaxCriterion::unpackElementData(CommunicationBuffer & buffer,
                                            const Array<Element> & elements,
                                            SynchronizationTag tag) {
  if (tag == SynchronizationTag::_nh_criterion) {
    this->unpackElementalDataHelper(criterion, buffer, elements, true,
                                    this->model.getFEEngine());
  }
}

/* -------------------------------------------------------------------------- */

} // namespace akantu

#endif /* __AKANTU_NEIGHBORHOOD_MAX_CRITERION_INLINE_IMPL_HH__ */
