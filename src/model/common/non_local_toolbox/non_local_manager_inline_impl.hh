/**
 * @file   non_local_manager_inline_impl.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Mon Sep 11 2017
 *
 * @brief  inline implementation of non-local manager functions
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
#include "neighborhood_base.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_NON_LOCAL_MANAGER_INLINE_IMPL_HH__
#define __AKANTU_NON_LOCAL_MANAGER_INLINE_IMPL_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
inline void NonLocalManager::registerNeighborhood(const ID & neighborhood,
                                                  const ID & weight_func_id) {

  /// check if neighborhood has already been created
  auto it = neighborhoods.find(neighborhood);
  if (it == neighborhoods.end()) {
    this->createNeighborhood(weight_func_id, neighborhood);
  }
}

/* -------------------------------------------------------------------------- */
inline NonLocalNeighborhoodBase &
NonLocalManager::getNeighborhood(const ID & name) const {
  AKANTU_DEBUG_IN();

  auto it = neighborhoods.find(name);

  AKANTU_DEBUG_ASSERT(it != neighborhoods.end(),
                      "The neighborhood " << name << " is not registered");

  AKANTU_DEBUG_OUT();
  return *(it->second);
}

} // namespace akantu

#endif /* __AKANTU_NON_LOCAL_MANAGER_INLINE_IMPL_HH__ */
