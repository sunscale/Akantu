/**
 * @file   neighborhood_base_inline_impl.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Feb 21 2013
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  Inline implementation of neighborhood base functions
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
#include "aka_grid_dynamic.hh"
#include "neighborhood_base.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_NEIGHBORHOOD_BASE_INLINE_IMPL_HH__
#define __AKANTU_NEIGHBORHOOD_BASE_INLINE_IMPL_HH__

namespace akantu {

inline void
NeighborhoodBase::insertIntegrationPoint(const IntegrationPoint & quad,
                                         const Vector<Real> & coords) {
  this->spatial_grid->insert(quad, coords);
}

} // namespace akantu

#endif /* __AKANTU_NEIGHBORHOOD_BASE_INLINE_IMPL_HH__ */
