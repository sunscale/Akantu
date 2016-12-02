/**
 * @file   node_synchronizer.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Sep 23 12:01:24 2016
 *
 * @brief  Implementation of the node synchronizer
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
#include "node_synchronizer.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
NodeSynchronizer::NodeSynchronizer(
    Mesh & mesh, const ID & id, MemoryID memory_id,
    const bool register_to_event_manager, StaticCommunicator & comm)
    : SynchronizerImpl<UInt>(id, memory_id, comm), mesh(mesh) {
  AKANTU_DEBUG_IN();

  if (register_to_event_manager) {
    this->mesh.registerEventHandler(*this);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
NodeSynchronizer::~NodeSynchronizer() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}


}  // akantu
