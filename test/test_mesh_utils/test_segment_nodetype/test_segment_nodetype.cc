/**
 * @file   test_segment_nodetype.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Sep 18 2015
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Test to verify that the node type is correctly associated to
 * the segments in parallel
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
#include "element_synchronizer.hh"
#include "mesh_utils.hh"

/* -------------------------------------------------------------------------- */
using namespace akantu;

int main(int argc, char * argv[]) {
  initialize(argc, argv);

  UInt spatial_dimension = 3;
  Mesh mesh(spatial_dimension);

  const auto & comm = Communicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  // partition the mesh
  if (prank == 0) {
    mesh.read("mesh.msh");
  }

  mesh.distribute();

  // compute the node types for each segment
  Mesh & mesh_facets = mesh.initMeshFacets();
  MeshUtils::buildSegmentToNodeType(mesh, mesh_facets);

  // verify that the number of segments per node type makes sense
  std::map<Int, UInt> nb_facets_per_nodetype;
  UInt nb_segments = 0;

  for (auto ghost_type : ghost_types) {
    const Array<Int> & segment_to_nodetype =
        mesh_facets.getData<Int>("segment_to_nodetype", _segment_2, ghost_type);

    // count the number of segments per node type
    for (auto & stn : segment_to_nodetype) {
      if (nb_facets_per_nodetype.find(stn) == nb_facets_per_nodetype.end())
        nb_facets_per_nodetype[stn] = 1;
      else
        ++nb_facets_per_nodetype[stn];
    }
    nb_segments += segment_to_nodetype.size();
  }

  // checking the solution
  if (nb_segments != 24)
    AKANTU_ERROR("The number of segments is wrong");

  if (prank == 0) {
    if (nb_facets_per_nodetype[-1] != 3 || nb_facets_per_nodetype[-2] != 9 ||
        nb_facets_per_nodetype[-3] != 12)
      AKANTU_ERROR("The segments of processor 0 have the wrong node type");

    if (nb_facets_per_nodetype.size() > 3)
      AKANTU_ERROR("Processor 0 cannot have any slave segment");
  }

  if (prank == psize - 1 &&
      nb_facets_per_nodetype.find(-2) != nb_facets_per_nodetype.end())
    AKANTU_ERROR("The last processor must not have any master facets");

  finalize();
  return 0;
}
