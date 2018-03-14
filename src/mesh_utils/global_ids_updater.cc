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
#include "global_ids_updater.hh"
#include "element_synchronizer.hh"
#include "mesh_accessor.hh"
#include "mesh_utils.hh"
/* -------------------------------------------------------------------------- */
#include <numeric>
/* -------------------------------------------------------------------------- */

namespace akantu {

UInt GlobalIdsUpdater::updateGlobalIDs(UInt local_nb_new_nodes) {
  UInt total_nb_new_nodes = this->updateGlobalIDsLocally(local_nb_new_nodes);

  if (mesh.isDistributed()) {
    this->synchronizeGlobalIDs();
  }
  return total_nb_new_nodes;
}

UInt GlobalIdsUpdater::updateGlobalIDsLocally(UInt local_nb_new_nodes) {
  const auto & comm = mesh.getCommunicator();
  Int rank = comm.whoAmI();
  Int nb_proc = comm.getNbProc();
  if (nb_proc == 1)
    return local_nb_new_nodes;

  /// resize global ids array
  Array<UInt> & nodes_global_ids = mesh.getGlobalNodesIds();
  UInt old_nb_nodes = mesh.getNbNodes() - local_nb_new_nodes;

  nodes_global_ids.resize(mesh.getNbNodes(), -1);

  /// compute the number of global nodes based on the number of old nodes
  Matrix<UInt> local_master_nodes(2, nb_proc, 0);
  for (UInt n = 0; n < old_nb_nodes; ++n)
    if (mesh.isLocalOrMasterNode(n))
      ++local_master_nodes(0, rank);

  /// compute amount of local or master doubled nodes
  for (UInt n = old_nb_nodes; n < mesh.getNbNodes(); ++n)
    if (mesh.isLocalOrMasterNode(n))
      ++local_master_nodes(1, rank);

  comm.allGather(local_master_nodes);

  local_master_nodes = local_master_nodes.transpose();
  UInt old_global_nodes =
      std::accumulate(local_master_nodes(0).storage(),
                      local_master_nodes(0).storage() + nb_proc, 0);

  /// update global number of nodes
  UInt total_nb_new_nodes =
      std::accumulate(local_master_nodes(1).storage(),
                      local_master_nodes(1).storage() + nb_proc, 0);

  if (total_nb_new_nodes == 0)
    return 0;

  /// set global ids of local and master nodes
  UInt starting_index =
      std::accumulate(local_master_nodes(1).storage(),
                      local_master_nodes(1).storage() + rank, old_global_nodes);

  for (UInt n = old_nb_nodes; n < mesh.getNbNodes(); ++n) {
    if (mesh.isLocalOrMasterNode(n)) {
      nodes_global_ids(n) = starting_index;
      ++starting_index;
    } else {
      nodes_global_ids(n) = -1;
    }
  }

  MeshAccessor mesh_accessor(mesh);
  mesh_accessor.setNbGlobalNodes(old_global_nodes + total_nb_new_nodes);
  return total_nb_new_nodes;
}

void GlobalIdsUpdater::synchronizeGlobalIDs() {
  this->reduce = true;
  this->synchronizer.slaveReductionOnce(*this, _gst_giu_global_conn);

#ifndef AKANTU_NDEBUG
  for (auto node : nodes_types) {
    auto node_type = mesh.getNodeType(node.first);
    if (node_type != _nt_pure_ghost)
      continue;
    auto n = 0u;

    for (auto & pair : node.second) {
      if (std::get<1>(pair) == _nt_pure_ghost)
        ++n;
    }

    if (n == node.second.size()) {
      AKANTU_DEBUG_WARNING(
          "The node " << n << "is ghost on all the neighboring processors");
    }
  }
#endif

  this->reduce = false;
  this->synchronizer.synchronizeOnce(*this, _gst_giu_global_conn);
}

} // akantu
