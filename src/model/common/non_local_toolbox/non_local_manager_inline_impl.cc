/**
 * @file   non_local_manager_inline_impl.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Mon Sep 21 17:30:03 2015
 *
 * @brief  inline implementation of non-local manager functions
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
__END_AKANTU__

#include "grid_synchronizer.hh"
#include "synchronizer_registry.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
inline void NonLocalManager::insertQuad(const QuadraturePoint & quad, const ID & id) {

  AKANTU_DEBUG_IN();
  const Vector<Real> coords = quad.getPosition();

  NeighborhoodMap::const_iterator it = neighborhoods.find(id);

  it->second->insertQuad(quad, coords);

  AKANTU_DEBUG_OUT();

}

/* -------------------------------------------------------------------------- */
inline NonLocalNeighborhoodBase & NonLocalManager::getNeighborhood(const ID & name) const{
  AKANTU_DEBUG_IN();
  ID tmp_name = name;
  if (name == "") tmp_name = default_neighborhood;

  NeighborhoodMap::const_iterator it = neighborhoods.find(tmp_name);

  AKANTU_DEBUG_ASSERT(it != neighborhoods.end(),
		      "The neighborhood " << tmp_name << " is not registered");

  AKANTU_DEBUG_OUT();
  return *(it->second);
}

/* -------------------------------------------------------------------------- */
inline void NonLocalManager::computeWeights() {
  AKANTU_DEBUG_IN();

  NeighborhoodMap::iterator it = neighborhoods.begin();
  NeighborhoodMap::iterator end = neighborhoods.end();

  for (; it != end; ++it)
    it->second->computeWeights();

  AKANTU_DEBUG_OUT();

}

/* -------------------------------------------------------------------------- */
inline void NonLocalManager::updatePairLists() {
  AKANTU_DEBUG_IN();

  NeighborhoodMap::iterator it = neighborhoods.begin();
  NeighborhoodMap::iterator end = neighborhoods.end();

  for (; it != end; ++it)
    it->second->updatePairList();

  AKANTU_DEBUG_OUT();

}








