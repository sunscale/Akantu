/**
 * @file   test_cohesive_parallel_intrinsic_implicit_insertion.cc
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 *
 * @date creation: Wed Nov 05 2014
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  Verifying the proper insertion and synchronization of intrinsic
 * cohesive elements
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include <fstream>
#include <iostream>
#include <limits>
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "material.hh"
#include "material_cohesive.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "mesh_utils.hh"
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */
using namespace akantu;
std::ofstream output;

/* -------------------------------------------------------------------------- */
void printMeshContent(Mesh & mesh) {

  const auto & comm = Communicator::getStaticCommunicator();
  Int prank = comm.whoAmI();
  comm.barrier();
  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {

    Mesh::type_iterator first =
        mesh.firstType(_all_dimensions, *gt, _ek_not_defined);
    Mesh::type_iterator last =
        mesh.lastType(_all_dimensions, *gt, _ek_not_defined);

    for (; first != last; ++first) {
      UInt nb_element = mesh.getNbElement(*first, *gt);
      output << std::endl
             << "Element type: " << *first << ", " << *gt << ": " << nb_element
             << " in the mesh of processor " << prank << std::endl;
      Array<UInt> & conn = mesh.getConnectivity(*first, *gt);
      for (UInt i = 0; i < conn.getSize(); ++i) {
        output << "Element no " << i << " ";
        for (UInt j = 0; j < conn.getNbComponent(); ++j) {
          output << conn(i, j) << " ";
        }
        output << std::endl;
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void printNodeList(Mesh & mesh) {

  Array<double> & nodes = mesh.getNodes();

  output << "Number of nodes: " << mesh.getNbNodes() << std::endl;
  for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
    output << "Node # " << i << ", x-coord: " << nodes(i, 0)
           << ", y-coord: " << nodes(i, 1)
           << ", of type: " << mesh.getNodeType(i) << std::endl;
  }
  output << std::endl;
}

/* -------------------------------------------------------------------------- */
void getGlobalIDs(Mesh & mesh) {

  const Array<UInt> & glob_id = mesh.getGlobalNodesIds();
  if (&glob_id) {
    output << "Global nodes ID: " << std::endl;
    for (UInt i = 0; i < glob_id.getSize(); ++i) {
      output << i << " " << glob_id(i) << std::endl;
    }
  }
  output << std::endl;
}

/* -------------------------------------------------------------------------- */
void printSynchroinfo(Mesh & mesh, const DistributedSynchronizer & synch) {

  const auto & comm = Communicator::getStaticCommunicator();
  Int prank = comm.whoAmI();
  Int psize = comm.getNbProc();

  if (comm.getNbProc() == 1)
    return;

  for (Int p = 0; p < psize; ++p) {

    if (p == prank)
      continue;
    output << "From processor " << prank << " to processor " << p << std::endl;
    const Array<Element> & sele = *(synch.getSendElement() + p);
    output << " Sending element(s): " << std::endl;
    for (UInt i = 0; i < sele.getSize(); ++i) {
      output << sele(i) << std::endl;
    }

    const Array<Element> & rele = *(synch.getReceiveElement() + p);
    output << " Receiving element(s): " << std::endl;
    for (UInt i = 0; i < rele.getSize(); ++i) {
      output << rele(i) << std::endl;
    }
  }
  output << std::endl;
}

/* -------------------------------------------------------------------------- */
void printDOF(SolidMechanicsModelCohesive & model) {

  const auto & comm = Communicator::getStaticCommunicator();
  if (comm.getNbProc() == 1)
    return;
  Int prank = comm.whoAmI();

  const DOFSynchronizer & dof = model.getDOFSynchronizer();

  output << "Number of global dofs " << dof.getNbGlobalDOFs()
         << " for processor " << prank << std::endl;

  const Array<UInt> & dof_global_ids = dof.getDOFGlobalIDs();

  for (UInt i = 0; i < dof_global_ids.getSize(); ++i) {

    output << "Local dof " << i << ", global id: " << dof_global_ids(i)
           << std::endl;
  }
  output << std::endl;
}

/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {

  std::string input_file = "input_file_iii.dat";
  std::string mesh_file = "2d_basic_interface.msh";
  std::string dir = "output_dir/";

  initialize(input_file, argc, argv);

  debug::setDebugLevel(dbl0);

  const UInt spatial_dimension = 2;

  Mesh mesh(spatial_dimension);

  const auto & comm = Communicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();
  akantu::MeshPartition * partition = NULL;

  std::stringstream filename;
  filename << dir.c_str() << "output_from_proc_" << prank << "_out_of_" << psize
           << ".out";
  output.open(filename.str());

  if (prank == 0) {

    mesh.read(mesh_file);
    partition = new MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);

    const ElementTypeMapArray<UInt> & partitions = partition->getPartitions();

    output << "The root processor read the mesh." << std::endl
           << "Only GMSH physical objects are created in the mesh."
           << std::endl;

    for (ghost_type_t::iterator gt = ghost_type_t::begin();
         gt != ghost_type_t::end(); ++gt) {

      Mesh::type_iterator first = mesh.firstType(_all_dimensions, *gt);
      Mesh::type_iterator last = mesh.lastType(_all_dimensions, *gt);

      for (; first != last; ++first) {

        output << "Element type: " << *first << " ghost type: " << *gt
               << std::endl;
        UInt nb_element = mesh.getNbElement(*first, *gt);
        output << nb_element << " to partitionate between " << psize
               << " processsors" << std::endl;

        Array<UInt> part = partitions(*first, *gt);

        for (UInt i = 0; i < part.getSize(); ++i) {
          output << i << " " << part(i) << std::endl;
        }
      }
    }

    output << "Nodes are also read and set with type -1 (normal node)"
           << std::endl;
    printNodeList(mesh);
  }

  SolidMechanicsModelCohesive model(mesh);

  output << "Before initParallel(), non-root processors have empty Mesh object"
         << std::endl;
  printMeshContent(mesh);

  model.initParallel(partition);

  output << "After initParallel(), Mesh object on each processor is a local "
            "partionated mesh containing ghost elements"
         << std::endl;
  printMeshContent(mesh);

  output << "Nodes are also partionated and new node types are defined:"
         << std::endl;
  printNodeList(mesh);
  output << "-3: pure ghost node -> not a local node" << std::endl
         << "-2: master node -> node shared with other processor(s) -> local "
            "and global node"
         << std::endl
         << ">0: slave node -> -> node shared with other processor(s) -> only "
            "local node (its id is the rank of the processor owning the master "
            "node)"
         << std::endl;

  output
      << "Each local node has a corresponding global id used during assembly: "
      << std::endl;
  getGlobalIDs(mesh);

  Mesh & mesh_facets = mesh.getMeshFacets();

  output << "Within cohesive element model, initParallel() creates a second "
            "Mesh object usually called mesh_facet"
         << std::endl
         << "This Mesh object contains all sub-dimensional elements where "
            "potential cohesive element can be inserted"
         << std::endl;
  printMeshContent(mesh_facets);

  const DistributedSynchronizer & synch_model = model.getSynchronizer();
  output << "The distributed synchronizer of solid mechanics model is used to "
            "synchronize fields with ghost element:"
         << std::endl;
  printSynchroinfo(mesh, synch_model);

  mesh.createGroupsFromMeshData<std::string>("physical_names");
  model.initFull(SolidMechanicsModelCohesiveOptions(_static));

  output << "In case of insertion along physical objects, cohesive elements "
            "are created during initFull()"
         << std::endl;
  output << "Elements list after insertion" << std::endl;
  printMeshContent(mesh);

  output << "Node list after insertion: (Total number of nodes "
         << mesh.getNbNodes() << ")" << std::endl;
  printNodeList(mesh);

  output << "Node global ids after insertion: (Total number of nodes "
         << mesh.getNbGlobalNodes() << ")" << std::endl;
  getGlobalIDs(mesh);

  const DistributedSynchronizer & coh_synch_model =
      *(model.getCohesiveSynchronizer());
  output << "Solid mechanics model cohesive has its own distributed "
            "synchronizer to handle ghost cohesive element:"
         << std::endl;
  printSynchroinfo(mesh, coh_synch_model);

  output << "A synchronizer dedicated to degrees of freedom (DOFs) is used by "
            "the solver to build matrices in parallel:"
         << std::endl
         << "This DOFSynchronizer is built based on nodes global id "
         << std::endl;
  printDOF(model);

  output.close();

  finalize();

  return EXIT_SUCCESS;
}
