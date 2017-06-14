/**
 * @file   test_segment_nodetype.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Sep 18 2015
 *
 * @brief  Test to verify that the node type is correctly associated to
 * the segments in parallel
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
#include "mesh_utils.hh"
#include "distributed_synchronizer.hh"

/* -------------------------------------------------------------------------- */
using namespace akantu;

int main(int argc, char * argv[]) {
  initialize(argc, argv);

  UInt spatial_dimension = 3;
  Mesh mesh(spatial_dimension);

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();
  DistributedSynchronizer * dist = NULL;

  // partition the mesh
  if (prank == 0) {
    mesh.read("mesh.msh");
    MeshPartitionScotch partition(mesh, spatial_dimension);
    partition.partitionate(psize);
    dist = DistributedSynchronizer::createDistributedSynchronizerMesh(
        mesh, &partition);
  } else {
    dist =
        DistributedSynchronizer::createDistributedSynchronizerMesh(mesh, NULL);
  }

  // compute the node types for each segment
  Mesh mesh_facets(mesh.initMeshFacets());
  MeshUtils::buildSegmentToNodeType(mesh, mesh_facets, dist);

  // verify that the number of segments per node type makes sense
  std::map<Int, UInt> nb_facets_per_nodetype;
  UInt nb_segments = 0;

  for (auto ghost_type : ghost_types) {
    const Array<Int> & segment_to_nodetype =
        mesh_facets.getData<Int>("segment_to_nodetype", _segment_2, ghost_type);

    // count the number of segments per node type
    for (UInt s = 0; s < segment_to_nodetype.getSize(); ++s) {
      if (nb_facets_per_nodetype.find(segment_to_nodetype(s)) ==
          nb_facets_per_nodetype.end())
        nb_facets_per_nodetype[segment_to_nodetype(s)] = 1;
      else
        ++nb_facets_per_nodetype[segment_to_nodetype(s)];
    }
    nb_segments += segment_to_nodetype.getSize();
  }

  // checking the solution
  if (nb_segments != 24)
    AKANTU_DEBUG_ERROR("The number of segments is wrong");

  if (prank == 0) {
    if (nb_facets_per_nodetype[-1] != 3 || nb_facets_per_nodetype[-2] != 9 ||
        nb_facets_per_nodetype[-3] != 12)
      AKANTU_DEBUG_ERROR(
          "The segments of processor 0 have the wrong node type");

    if (nb_facets_per_nodetype.size() > 3)
      AKANTU_DEBUG_ERROR("Processor 0 cannot have any slave segment");
  }

  if (prank == psize - 1 &&
      nb_facets_per_nodetype.find(-2) != nb_facets_per_nodetype.end())
    AKANTU_DEBUG_ERROR("The last processor must not have any master facets");

  finalize();
  return 0;
}
