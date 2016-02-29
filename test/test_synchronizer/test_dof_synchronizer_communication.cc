/**
 * @file   test_dof_synchronizer_communication.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 * @date creation: Tue Dec 09 2014
 *
 * @brief  test to synchronize global equation numbers
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "distributed_synchronizer.hh"
#include "dof_synchronizer.hh"
#include "synchronizer_registry.hh"
#include "mesh.hh"
#include "mesh_partition_scotch.hh"
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#  include "dumper_paraview.hh"
#endif //AKANTU_USE_IOHELPER

#include "test_dof_data_accessor.hh"


using namespace akantu;

/* -------------------------------------------------------------------------- */
/* Main                                                                       */
/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
  initialize(argc, argv);

  UInt spatial_dimension = 3;

  Mesh mesh(spatial_dimension);

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();
  bool wait=true;
  if(argc >1) {
    if(prank == 0)
    while(wait);
  }

  DistributedSynchronizer * communicator = NULL;
  if(prank == 0) {
    mesh.read("cube.msh");

    MeshPartition * partition = new MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);
    communicator = DistributedSynchronizer::createDistributedSynchronizerMesh(mesh, partition);
    delete partition;
  } else {
    communicator = DistributedSynchronizer::createDistributedSynchronizerMesh(mesh, NULL);
  }

/* -------------------------------------------------------------------------- */
/* test the communications of the dof synchronizer                                */
/* -------------------------------------------------------------------------- */

  std::cout << "Initializing the synchronizer" << std::endl;
  DOFSynchronizer dof_synchronizer(mesh, spatial_dimension);
  dof_synchronizer.initGlobalDOFEquationNumbers();

  AKANTU_DEBUG_INFO("Creating TestDOFAccessor");
  TestDOFAccessor test_dof_accessor(dof_synchronizer.getGlobalDOFEquationNumbers());
  SynchronizerRegistry synch_registry(test_dof_accessor);
  synch_registry.registerSynchronizer(dof_synchronizer,_gst_test);

  AKANTU_DEBUG_INFO("Synchronizing tag");
  synch_registry.synchronize(_gst_test);

  delete communicator;
  finalize();

  return EXIT_SUCCESS;
}

