/**
 * @file   periodic_node_synchronizer.cc
 *
 * @author Nicolas Richart
 *
 * @date creation  Tue May 29 2018
 *
 * @brief Implementation of the periodic node synchronizer
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
#include "periodic_node_synchronizer.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
PeriodicNodeSynchronizer::PeriodicNodeSynchronizer(
    Mesh & mesh, const ID & id, MemoryID memory_id,
    const bool register_to_event_manager, EventHandlerPriority event_priority)
    : NodeSynchronizer(mesh, id + ":masters", memory_id,
                       register_to_event_manager, event_priority) {}

/* -------------------------------------------------------------------------- */
void PeriodicNodeSynchronizer::update() {
  this->copySchemes(this->mesh.getNodeSynchronizer());

  this->filterScheme([&](auto && node) { return mesh.isPeriodicMaster(node); });
}

/* -------------------------------------------------------------------------- */
void PeriodicNodeSynchronizer::synchronizeOnceImpl(
    DataAccessor<UInt> & data_accessor, const SynchronizationTag & tag) const {}

/* -------------------------------------------------------------------------- */
void PeriodicNodeSynchronizer::asynchronousSynchronizeImpl(
    const DataAccessor<UInt> & data_accessor,
    const SynchronizationTag & tag) {}

/* -------------------------------------------------------------------------- */
void PeriodicNodeSynchronizer::waitEndSynchronizeImpl(
    DataAccessor<UInt> & data_accessor,
    const SynchronizationTag & tag) {}

} // namespace akantu
