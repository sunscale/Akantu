/**
 * @file   global_ids_updater.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Apr 13 2012
 * @date last modification: Fri Dec 08 2017
 *
 * @brief  Functions of the GlobalIdsUpdater
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
#include "global_ids_updater.hh"
#include "element_synchronizer.hh"
#include "mesh_accessor.hh"
#include "mesh_utils.hh"
/* -------------------------------------------------------------------------- */
#include <numeric>
/* -------------------------------------------------------------------------- */

namespace akantu {

UInt GlobalIdsUpdater::updateGlobalIDs(UInt local_nb_new_nodes) {
  if (mesh.getCommunicator().getNbProc() == 1) {
    return local_nb_new_nodes;
  }

  UInt total_nb_new_nodes = this->updateGlobalIDsLocally(local_nb_new_nodes);

  if (mesh.isDistributed()) {
    this->synchronizeGlobalIDs();
  }
  return total_nb_new_nodes;
}

UInt GlobalIdsUpdater::updateGlobalIDsLocally(UInt local_nb_new_nodes) {
  const auto & comm = mesh.getCommunicator();
  Int nb_proc = comm.getNbProc();
  if (nb_proc == 1) {
    return local_nb_new_nodes;
  }

  /// resize global ids array
  MeshAccessor mesh_accessor(mesh);
  auto && nodes_global_ids = mesh_accessor.getNodesGlobalIds();
  UInt old_nb_nodes = mesh.getNbNodes() - local_nb_new_nodes;

  nodes_global_ids.resize(mesh.getNbNodes(), -1);

  /// compute the number of global nodes based on the number of old nodes
  Vector<UInt> local_master_nodes(2, 0);
  auto range_old = arange(old_nb_nodes);
  local_master_nodes(0) =
      aka::count_if(range_old.begin(), range_old.end(),
                    [&](auto && n) { return mesh.isLocalOrMasterNode(n); });

  /// compute amount of local or master doubled nodes
  auto range_new = arange(old_nb_nodes, mesh.getNbNodes());
  local_master_nodes(1) =
      aka::count_if(range_new.begin(), range_new.end(),
                    [&](auto && n) { return mesh.isLocalOrMasterNode(n); });

  auto starting_index = local_master_nodes(1);

  comm.allReduce(local_master_nodes);

  UInt old_global_nodes = local_master_nodes(0);
  UInt total_nb_new_nodes = local_master_nodes(1);

  if (total_nb_new_nodes == 0) {
    return 0;
  }

  /// set global ids of local and master nodes
  comm.exclusiveScan(starting_index);
  starting_index += old_global_nodes;

  for (auto n : range_new) {
    if (mesh.isLocalOrMasterNode(n)) {
      nodes_global_ids(n) = starting_index;
      ++starting_index;
    }
  }

  mesh_accessor.setNbGlobalNodes(old_global_nodes + total_nb_new_nodes);
  return total_nb_new_nodes;
}

void GlobalIdsUpdater::synchronizeGlobalIDs() {
  this->reduce = true;
  this->synchronizer.slaveReductionOnce(*this,
                                        SynchronizationTag::_giu_global_conn);

#ifndef AKANTU_NDEBUG
  for (auto node : nodes_flags) {
    auto node_flag = mesh.getNodeFlag(node.first);
    if (node_flag != NodeFlag::_pure_ghost) {
      continue;
    }
    auto n = 0U;

    for (auto & pair : node.second) {
      if (std::get<1>(pair) == NodeFlag::_pure_ghost) {
        ++n;
      }
    }

    if (n == node.second.size()) {
      AKANTU_DEBUG_WARNING(
          "The node " << n << "is ghost on all the neighboring processors");
    }
  }
#endif

  this->reduce = false;
  this->synchronizer.synchronizeOnce(*this,
                                     SynchronizationTag::_giu_global_conn);
}

} // namespace akantu
