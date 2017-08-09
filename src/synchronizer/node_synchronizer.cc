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
NodeSynchronizer::NodeSynchronizer(Mesh & mesh, const ID & id,
                                   MemoryID memory_id,
                                   const bool register_to_event_manager,
                                   EventHandlerPriority event_priority)
    : SynchronizerImpl<UInt>(id, memory_id, mesh.getCommunicator()),
      mesh(mesh) {
  AKANTU_DEBUG_IN();

  if (register_to_event_manager) {
    this->mesh.registerEventHandler(*this, event_priority);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
NodeSynchronizer::~NodeSynchronizer() = default;

/* -------------------------------------------------------------------------- */
void NodeSynchronizer::onNodesAdded(const Array<UInt> & nodes_list,
                                    const NewNodesEvent &) {
  std::map<UInt, std::vector<UInt>> nodes_per_proc;
  Array<UInt> sizes_per_proc(nb_proc, nb_proc, UInt(0));
  Vector<UInt> local_sizes_per_proc = sizes_per_proc.begin(nb_proc)[rank];

  for (auto & local_id : nodes_list) {
    auto type = mesh.getNodeType(local_id);
    if (type <= 0)
      continue; // local, master or pure ghost

    auto global_id = mesh.getNodeGlobalId(local_id);
    auto proc = UInt(type);
    nodes_per_proc[proc].push_back(global_id);
    ++local_sizes_per_proc[proc];

    auto & scheme = communications.getScheme(proc, _recv);
    scheme.push_back(local_id);
  }

  communicator.allGather(sizes_per_proc);
  std::vector<CommunicationRequest> send_request_per_proc,
      recv_request_per_proc;
  std::map<UInt, std::vector<UInt>> nodes_needed_by_proc;
  for (UInt proc = 0; proc < nb_proc; ++proc) {
    auto size = sizes_per_proc(proc, rank);
    if (size == 0)
      continue;

    nodes_needed_by_proc[proc].resize(size);
    recv_request_per_proc.push_back(communicator.asyncReceive(
        nodes_needed_by_proc[proc], proc, Tag::genTag(rank, 0, 0)));
  }


  for (auto & pair : nodes_per_proc) {
    auto proc = pair.first;
    auto & nodes = pair.second;

    send_request_per_proc.push_back(communicator.asyncSend(
        nodes, proc, Tag::genTag(proc, 0, 0)));
  }

  UInt req_nb;
  while((req_nb = communicator.waitAny(recv_request_per_proc)) != UInt(-1)) {
    auto & request = recv_request_per_proc[req_nb];
    auto proc = request.getSource();

    auto & nodes = nodes_needed_by_proc[proc];

    auto & scheme = communications.getScheme(proc, _send);
    for (auto global_id : nodes) {
      auto local_id = mesh.getNodeLocalId(global_id);
      scheme.push_back(local_id);
    }
  }

  communicator.waitAll(send_request_per_proc);
  communicator.freeCommunicationRequest(send_request_per_proc);
}

} // namespace akantu
