/**
 * @file   test_dof_synchronizer.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 17 2011
 * @date last modification: Thu Apr 03 2014
 *
 * @brief  Test the functionality of the DOFSynchronizer class
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_common.hh"
#include "dof_synchronizer.hh"
#include "mesh_partition_scotch.hh"
#include "mesh_io.hh"
#include "static_communicator.hh"
#include "distributed_synchronizer.hh"

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.hh"
#endif //AKANTU_USE_IOHELPER

/* -------------------------------------------------------------------------- */
using namespace akantu;

int main(int argc, char *argv[])
{
  const UInt spatial_dimension = 2;

  initialize(argc, argv);

  StaticCommunicator & comm = akantu::StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  Mesh mesh(spatial_dimension);

  /* ------------------------------------------------------------------------ */
  /* Parallel initialization                                                  */
  /* ------------------------------------------------------------------------ */
  DistributedSynchronizer * communicator;
  MeshPartition * partition;

  if(prank == 0) {
    MeshIOMSH mesh_io;
    mesh_io.read("bar.msh", mesh);

    std::cout << "Partitioning mesh..." << std::endl;
    partition = new akantu::MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);
    communicator = DistributedSynchronizer::createDistributedSynchronizerMesh(mesh, partition);
    delete partition;
  } else {
    communicator = DistributedSynchronizer::createDistributedSynchronizerMesh(mesh, NULL);
  }

  UInt nb_nodes = mesh.getNbNodes();

  Array<Real> dof_vector(nb_nodes, spatial_dimension, "Test vector");

  std::cout << "Initializing the synchronizer" << std::endl;
  DOFSynchronizer dof_synchronizer(mesh, spatial_dimension);


  /* ------------------------------------------------------------------------ */
  /* test the sznchroniyation                                                 */
  /* ------------------------------------------------------------------------ */
  for (UInt n = 0; n < nb_nodes; ++n) {
    UInt gn = mesh.getNodeGlobalId(n);
    for (UInt d = 0; d < spatial_dimension; ++d) {
      if(mesh.isMasterNode(n))     dof_vector(n,d) = gn*spatial_dimension + d;
      else if(mesh.isLocalNode(n)) dof_vector(n,d) = - (double) (gn*spatial_dimension + d);
      else if(mesh.isSlaveNode(n)) dof_vector(n,d) = NAN;
      else dof_vector(n,d) = -NAN;
    }
  }

  std::cout << "Synchronizing a dof vector" << std::endl;
  dof_synchronizer.synchronize(dof_vector);

  for (UInt n = 0; n < nb_nodes; ++n) {
    UInt gn = mesh.getNodeGlobalId(n);
    for (UInt d = 0; d < spatial_dimension; ++d) {
      if(!((mesh.isMasterNode(n) && dof_vector(n,d) == (gn*spatial_dimension + d)) ||
	   (mesh.isLocalNode(n) && dof_vector(n,d) == - (double) (gn*spatial_dimension + d)) ||
	   (mesh.isSlaveNode(n) && dof_vector(n,d) == (gn*spatial_dimension + d)) ||
	   (mesh.isPureGhostNode(n))
	   )
	 )
	{
	  debug::setDebugLevel(dblTest);
	  std::cout << "prank : " << prank << " (node " << gn*spatial_dimension + d << "[" << n * spatial_dimension + d << "]) - "
		    << (mesh.isMasterNode(n) && dof_vector(n,d) == (gn*spatial_dimension + d)) << " "
		    << (mesh.isLocalNode(n) && dof_vector(n,d) == - (double) (gn*spatial_dimension + d)) << " "
		    << (mesh.isSlaveNode(n) && dof_vector(n,d) == (gn*spatial_dimension + d)) << " "
		    << (mesh.isPureGhostNode(n)) << std::endl;
	  std::cout << dof_vector << dof_synchronizer.getDOFGlobalIDs() << dof_synchronizer.getDOFTypes() << std::endl;
	  debug::setDebugLevel(dblDebug);
	  return EXIT_FAILURE;
	}
    }
  }


  /* ------------------------------------------------------------------------ */
  /* test the reduce operation                                                */
  /* ------------------------------------------------------------------------ */
  for (UInt n = 0; n < nb_nodes; ++n) {
    for (UInt d = 0; d < spatial_dimension; ++d) {
      if(mesh.isMasterNode(n))     dof_vector(n,d) = 1;
      else if(mesh.isLocalNode(n)) dof_vector(n,d) = -300;
      else if(mesh.isSlaveNode(n)) dof_vector(n,d) = 2;
      else dof_vector(n,d) = -500;
    }
  }

  std::cout << "Reducing a dof vector" << std::endl;
  dof_synchronizer.reduceSynchronize<AddOperation>(dof_vector);

  for (UInt n = 0; n < nb_nodes; ++n) {
    for (UInt d = 0; d < spatial_dimension; ++d) {
      if(!((mesh.isMasterNode(n) && dof_vector(n,d) >= 3) ||
	   (mesh.isLocalNode(n) && dof_vector(n,d) == -300) ||
	   (mesh.isSlaveNode(n) && dof_vector(n,d) >= 3) ||
	   (mesh.isPureGhostNode(n) && dof_vector(n,d) == -500)
	   )
	 )
	{
	  debug::setDebugLevel(dblTest);
	  std::cout << dof_vector
		    << dof_synchronizer.getDOFGlobalIDs()
		    << dof_synchronizer.getDOFTypes() << std::endl;
	  debug::setDebugLevel(dblDebug);
	  return EXIT_FAILURE;
	}
    }
  }

  /* ------------------------------------------------------------------------ */
  /* test the gather/scatter                                                  */
  /* ------------------------------------------------------------------------ */
  dof_vector.clear();
  for (UInt n = 0; n < nb_nodes; ++n) {
    UInt gn = mesh.getNodeGlobalId(n);
    for (UInt d = 0; d < spatial_dimension; ++d) {
      if(mesh.isMasterNode(n))     dof_vector(n,d) = gn * spatial_dimension + d;
      else if(mesh.isLocalNode(n)) dof_vector(n,d) = - (double) (gn * spatial_dimension + d);
      else if(mesh.isSlaveNode(n)) dof_vector(n,d) = NAN;
      else dof_vector(n,d) = -NAN;
    }
  }

  std::cout << "Initializing the gather/scatter information" << std::endl;
  dof_synchronizer.initScatterGatherCommunicationScheme();

  std::cout << "Gathering on proc 0" << std::endl;
  if(prank == 0) {
    UInt nb_global_nodes = mesh.getNbGlobalNodes();
    Array<Real> gathered(nb_global_nodes, spatial_dimension, "gathered information");
    dof_synchronizer.gather(dof_vector, 0, &gathered);
    for (UInt n = 0; n < nb_nodes; ++n) {
      for (UInt d = 0; d < spatial_dimension; ++d) {
	if(std::abs(gathered(n,d)) != n * spatial_dimension + d) {
	  debug::setDebugLevel(dblTest);
	  std::cout << gathered << std::endl;
	  std::cout << dof_vector
		    << dof_synchronizer.getDOFGlobalIDs()
		    << dof_synchronizer.getDOFTypes() << std::endl;
	  debug::setDebugLevel(dblDebug);
	  return EXIT_FAILURE;
	}
      }
    }
  } else {
    dof_synchronizer.gather(dof_vector, 0);
  }

  dof_vector.clear();
  std::cout << "Scattering from proc 0" << std::endl;
  if(prank == 0) {
    UInt nb_global_nodes = mesh.getNbGlobalNodes();
    Array<Real> to_scatter(nb_global_nodes, spatial_dimension, "to scatter information");
    for (UInt d = 0; d < nb_global_nodes * spatial_dimension; ++d) {
      to_scatter.storage()[d] = d;
    }
    dof_synchronizer.scatter(dof_vector, 0, &to_scatter);
  } else {
    dof_synchronizer.scatter(dof_vector, 0);
  }

  for (UInt n = 0; n < nb_nodes; ++n) {
    UInt gn = mesh.getNodeGlobalId(n);
    for (UInt d = 0; d < spatial_dimension; ++d) {
      if(!mesh.isPureGhostNode(n) && !(dof_vector(n,d) == (gn * spatial_dimension + d))) {
	debug::setDebugLevel(dblTest);
	std::cout << dof_vector
		  << dof_synchronizer.getDOFGlobalIDs()
		  << dof_synchronizer.getDOFTypes() << std::endl;
	debug::setDebugLevel(dblDebug);
	return EXIT_FAILURE;
      }
    }
  }

  delete communicator;
  finalize();

  return 0;
}
