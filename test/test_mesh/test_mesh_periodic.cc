/**
 * @file   test_mesh_periodic.cc
 *
 * @author Nicolas Richart
 *
 * @date creation  Sun Feb 11 2018
 *
 * @brief test makePeriodic
 *
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
#include "data_accessor.hh"
#include "mesh.hh"
#include "mesh_accessor.hh"
//#include "mesh_partition_scotch.hh"
#include "periodic_node_synchronizer.hh"
/* -------------------------------------------------------------------------- */
#include "dumpable_inline_impl.hh"
//#include "dumper_element_partition.hh"
#include "dumper_iohelper_paraview.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char ** argv) {
  initialize(argc, argv);

  constexpr UInt dim = 3;

  auto prank = Communicator::getStaticCommunicator().whoAmI();
  // auto psize = Communicator::getStaticCommunicator().getNbProc();

  Mesh mesh(dim);
  if (prank == 0) {
    mesh.read("cube_periodic.msh");
  }

  MeshAccessor mesh_accessor(mesh);
  // mesh_accessor.wipePeriodicInfo();
  // mesh.makePeriodic(_z);

  // if (prank == 0) {
  //   MeshPartitionScotch partition(mesh, dim);
  //   partition.partitionate(psize);
  // }

  UInt offset = 0;
  for (auto && type : mesh.elementTypes()) {
    auto & g_ids = mesh.getDataPointer<UInt>("global_ids", type);
    for (auto && data : enumerate(g_ids)) {
      std::get<1>(data) = offset + std::get<0>(data);
    }
    offset += g_ids.size();
  }

  mesh.distribute();

  mesh.makePeriodic(_x);
  mesh.makePeriodic(_y);
  mesh.makePeriodic(_z);

  auto * dumper = new DumperParaview("periodic", "./paraview");
  mesh.registerExternalDumper(*dumper, "periodic", true);
  mesh.addDumpMesh(mesh);
  if (mesh.isDistributed()) {
    mesh.addDumpFieldExternalToDumper(
        "periodic", "node_type",
        const_cast<const Mesh &>(mesh).getNodesFlags());
  }
  mesh.dump();

  Array<Int> periodic(mesh.getNbNodes(), 1, 0.);
  Array<Int> masters(mesh.getNbNodes(), 1, 0.);
  Array<Int> global_ids(mesh.getNbNodes(), 1, 0.);
  UInt prev_node = -1;
  UInt value = 0;
  const auto & periodic_ms = mesh.getPeriodicMasterSlaves();
  for (auto & pair : periodic_ms) {
    if (prev_node != pair.first) {
      ++value;
    }

    prev_node = pair.first;
    periodic(pair.first) = value;
    periodic(pair.second) = value;

    masters(pair.first) = 1;
    global_ids(pair.first) = mesh.getNodeGlobalId(pair.second);

    auto it = periodic_ms.find(pair.second);
    if (it != periodic_ms.end()) {
      AKANTU_EXCEPTION(pair.second << " is slave of " << pair.first
                                   << " and master of " << it->second);
    }
  }
  mesh.addDumpFieldExternalToDumper("periodic", "periodic", periodic);
  mesh.addDumpFieldExternalToDumper("periodic", "masters", masters);
  mesh.addDumpFieldExternalToDumper("periodic", "global_ids", global_ids);
  mesh.addDumpFieldExternalToDumper("periodic", "element_global_ids",
                                    mesh.getData<UInt>("global_ids"));

  mesh.dump();

  Array<Int> data(mesh.getNbNodes(), 1, 0.);
  mesh.addDumpFieldExternalToDumper("periodic", "data", data);
  for (auto node : arange(mesh.getNbNodes())) {
    if (mesh.isPeriodicMaster(node)) {
      data(node) = 1 * (prank + 1);
      if (mesh.isMasterNode(node) or mesh.isLocalNode(node)) {
        data(node) = 10 * (prank + 1);
      }
    }
  }

  mesh.dump();

  // SimpleUIntDataAccessor<Int> data_accessor(data,
  // SynchronizationTag::_user_1);
  // mesh.getPeriodicNodeSynchronizer().synchronizeOnce(data_accessor,
  //                                                    SynchronizationTag::_user_1);
  mesh.dump();
}
