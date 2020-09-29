/**
 * @file   mesh_utils_distribution.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Nov 08 2016
 * @date last modification: Tue Nov 07 2017
 *
 * @brief  Implementation of the methods of mesh  utils distribute
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
#include "mesh_utils_distribution.hh"
#include "element_info_per_processor.hh"
#include "element_synchronizer.hh"
#include "mesh.hh"
#include "mesh_accessor.hh"
#include "mesh_partition.hh"
#include "mesh_utils.hh"
#include "node_info_per_processor.hh"
#include "node_synchronizer.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

namespace MeshUtilsDistribution {
  /* ------------------------------------------------------------------------ */
  void distributeMeshCentralized(Mesh & mesh, UInt /*unused*/,
                                 const MeshPartition & partition) {
    MeshAccessor mesh_accessor(mesh);
    ElementSynchronizer & element_synchronizer =
        mesh_accessor.getElementSynchronizer();
    NodeSynchronizer & node_synchronizer = mesh_accessor.getNodeSynchronizer();

    const Communicator & comm = element_synchronizer.getCommunicator();

    UInt nb_proc = comm.getNbProc();
    UInt my_rank = comm.whoAmI();

    mesh_accessor.setNbGlobalNodes(mesh.getNbNodes());
    auto & gids = mesh_accessor.getNodesGlobalIds();

    if (nb_proc == 1) {
      return;
    }

    gids.resize(0);

    mesh.synchronizeGroupNames();

    AKANTU_DEBUG_ASSERT(
        partition.getNbPartition() == nb_proc,
        "The number of partition does not match the number of processors: "
            << partition.getNbPartition() << " != " << nb_proc);

    /**
     * connectivity and communications scheme construction
     */
    UInt count = 0;
    /* --- MAIN LOOP ON TYPES --- */
    for (auto && type :
         mesh.elementTypes(_all_dimensions, _not_ghost, _ek_not_defined)) {
      /// \todo change this ugly way to avoid a problem if an element
      /// type is present in the mesh but not in the partitions
      try {
        partition.getPartition(type, _not_ghost);
      } catch (...) {
        continue;
      }

      MasterElementInfoPerProc proc_infos(element_synchronizer, count, my_rank,
                                          type, partition);
      proc_infos.synchronize();
      ++count;
    }

    { /// Ending the synchronization of elements by sending a stop message
      MasterElementInfoPerProc proc_infos(element_synchronizer, count, my_rank,
                                          _not_defined, partition);
      proc_infos.synchronize();
      ++count;
    }

    /**
     * Nodes synchronization
     */
    MasterNodeInfoPerProc node_proc_infos(node_synchronizer, count, my_rank);
    node_proc_infos.synchronize();

    MeshUtils::fillElementToSubElementsData(mesh);

    mesh_accessor.setDistributed();

    AKANTU_DEBUG_OUT();
  }

  /* ------------------------------------------------------------------------ */
  void distributeMeshCentralized(Mesh & mesh, UInt root) {
    MeshAccessor mesh_accessor(mesh);
    ElementSynchronizer & element_synchronizer =
        mesh_accessor.getElementSynchronizer();
    NodeSynchronizer & node_synchronizer = mesh_accessor.getNodeSynchronizer();

    const Communicator & comm = element_synchronizer.getCommunicator();

    UInt nb_proc = comm.getNbProc();

    mesh_accessor.getNodesGlobalIds().resize(0);

    if (nb_proc == 1) {
      return;
    }

    mesh.synchronizeGroupNames();

    /**
     * connectivity and communications scheme construction on distant
     * processors
     */
    UInt count = 0;
    bool need_synchronize = true;
    do {
      /* --------<<<<-SIZE--------------------------------------------------- */
      SlaveElementInfoPerProc proc_infos(element_synchronizer, count, root);
      need_synchronize = proc_infos.synchronize();

      ++count;
    } while (need_synchronize);

    /**
     * Nodes synchronization
     */

    SlaveNodeInfoPerProc node_proc_infos(node_synchronizer, count, root);
    node_proc_infos.synchronize();

    MeshUtils::fillElementToSubElementsData(mesh);

    mesh_accessor.setDistributed();
  }

} // namespace MeshUtilsDistribution

} // namespace akantu
