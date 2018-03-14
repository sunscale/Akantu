/**
* @file   node_synchronizer.cc
*
* @author Nicolas Richart <nicolas.richart@epfl.ch>
*
* @date creation: Fri Jun 18 2010
* @date last modification: Wed Nov 15 2017
*
* @brief  Implementation of the node synchronizer
*
* @section LICENSE
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
#include "node_synchronizer.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
NodeSynchronizer::NodeSynchronizer(Mesh & mesh, const ID & id,
                                   MemoryID memory_id,
                                   const bool register_to_event_manager,
                                   EventHandlerPriority event_priority)
    : SynchronizerImpl<UInt>(mesh.getCommunicator(), id, memory_id),
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
void NodeSynchronizer::onNodesAdded(const Array<UInt> & /*nodes_list*/,
                                    const NewNodesEvent &) {
  std::map<UInt, std::vector<UInt>> nodes_per_proc;

  // recreates fully the schemes due to changes of global ids
  // \TODO add an event to handle global id changes
  for(auto && data : communications.iterateRecvSchemes()) {
    auto & scheme = data.second;
    scheme.resize(0);
  }

  for (auto && local_id : arange(mesh.getNbNodes())) {
    auto type = mesh.getNodeType(local_id);
    if (type < 0)
      continue; // local, master or pure ghost

    auto global_id = mesh.getNodeGlobalId(local_id);
    auto proc = UInt(type);
    nodes_per_proc[proc].push_back(global_id);

    auto & scheme = communications.createScheme(proc, _recv);
    scheme.push_back(local_id);
  }

  std::vector<CommunicationRequest> send_requests;
  for (auto & pair : nodes_per_proc) {
    auto proc = pair.first;
    auto & nodes = pair.second;

    send_requests.push_back(
        communicator.asyncSend(nodes, proc, Tag::genTag(proc, 0, 0xcafe)));
  }

  Array<UInt> buffer;

  communicator.receiveAnyNumber(
      send_requests, buffer,
      [&](auto && proc, auto && nodes) {
        auto & scheme = communications.createScheme(proc, _send);
        scheme.resize(nodes.size());
        for (auto && data : enumerate(nodes)) {
          auto global_id = std::get<1>(data);
          auto local_id = mesh.getNodeLocalId(global_id);
          AKANTU_DEBUG_ASSERT(local_id != UInt(-1),
                              "The global node " << global_id
                                                 << "is not known on rank "
                                                 << rank);
          scheme[std::get<0>(data)] = local_id;
        }
      },
      Tag::genTag(rank, 0, 0xcafe));
}

/* -------------------------------------------------------------------------- */
UInt NodeSynchronizer::sanityCheckDataSize(const Array<UInt> & nodes,
                                           const SynchronizationTag & tag,
                                           bool from_comm_desc) const {
  UInt size =
      SynchronizerImpl<UInt>::sanityCheckDataSize(nodes, tag, from_comm_desc);

  // positions
  size += mesh.getSpatialDimension() * sizeof(Real) * nodes.size();

  return size;
}

/* -------------------------------------------------------------------------- */
void NodeSynchronizer::packSanityCheckData(
    CommunicationBuffer & buffer, const Array<UInt> & nodes,
    const SynchronizationTag & /*tag*/) const {
  auto dim = mesh.getSpatialDimension();
  for (auto && node : nodes) {
    buffer << Vector<Real>(mesh.getNodes().begin(dim)[node]);
  }
}

/* -------------------------------------------------------------------------- */
void NodeSynchronizer::unpackSanityCheckData(CommunicationBuffer & buffer,
                                             const Array<UInt> & nodes,
                                             const SynchronizationTag & tag,
                                             UInt proc, UInt rank) const {
  auto dim = mesh.getSpatialDimension();
  // std::set<SynchronizationTag> skip_conn_tags{_gst_smmc_facets_conn,
  //                                             _gst_giu_global_conn};

  // bool is_skip_tag_conn = skip_conn_tags.find(tag) != skip_conn_tags.end();

  for (auto && node : nodes) {
    Vector<Real> pos_remote(dim);
    buffer >> pos_remote;
    Vector<Real> pos(mesh.getNodes().begin(dim)[node]);

    auto dist = pos_remote.distance(pos);
    if (not Math::are_float_equal(dist, 0.)) {
      AKANTU_EXCEPTION("Unpacking an unknown value for the node "
                       << node << "(position " << pos << " != buffer "
                       << pos_remote << ") [" << dist << "] - tag: " << tag
                       << " comm from " << proc << " to " << rank);
    }
  }
}
/* -------------------------------------------------------------------------- */

} // namespace akantu
