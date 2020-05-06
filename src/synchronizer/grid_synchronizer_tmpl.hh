/**
 * @file   grid_synchronizer_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Jul 06 2017
 * @date last modification: Wed Aug 09 2017
 *
 * @brief  implementation of the templated part of the grid syncrhonizers
 *
 *
 * Copyright (©) 2016-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "grid_synchronizer.hh"

#ifndef __AKANTU_GRID_SYNCHRONIZER_TMPL_HH__
#define __AKANTU_GRID_SYNCHRONIZER_TMPL_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
template <typename E>
GridSynchronizer::GridSynchronizer(Mesh & mesh, const SpatialGrid<E> & grid,
                                   const ID & id, MemoryID memory_id,
                                   const bool register_to_event_manager,
                                   EventHandlerPriority event_priority)
    : ElementSynchronizer(mesh, id, memory_id, register_to_event_manager,
                          event_priority) {
  AKANTU_DEBUG_IN();

  this->createGridSynchronizer(grid);

  AKANTU_DEBUG_OUT();
}

template <typename E>
GridSynchronizer::GridSynchronizer(
    Mesh & mesh, const SpatialGrid<E> & grid,
    SynchronizerRegistry & synchronizer_registry,
    const std::set<SynchronizationTag> & tags_to_register, const ID & id,
    MemoryID memory_id, const bool register_to_event_manager,
    EventHandlerPriority event_priority)
    : GridSynchronizer(mesh, grid, id, memory_id, register_to_event_manager,
                       event_priority) {
  AKANTU_DEBUG_IN();

  // Register the tags if any
  for (auto & tag : tags_to_register) {
    synchronizer_registry.registerSynchronizer(*this, tag);
  }

  AKANTU_DEBUG_OUT();
}

} // namespace akantu

#endif /* __AKANTU_GRID_SYNCHRONIZER_TMPL_HH__ */
