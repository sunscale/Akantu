/**
 * @file   test_mesh_boundary.cc
 *
 *
 * @date creation: Fri May 03 2013
 * @date last modification: Tue Sep 02 2014
 *
 * @brief  Thest the element groups 
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#include <iostream>
#include <sstream>
#include "aka_common.hh"
#include "mesh.hh"

using namespace akantu;

/* -------------------------------------------------------------------------- */

int main(int argc, char* argv[]) {
  UInt spatialDimension(3);

  akantu::initialize(argc, argv);

  Mesh mesh(spatialDimension, "mesh_names");

  std::cout << "Loading the mesh." << std::endl;

  //    mesh.read("./cube_physical_names.msh");
  mesh.read("./cube_physical_names.msh");
  std::stringstream sstr;

  std::cout << "Examining mesh:" << std::endl;

  // Inspection of the number of boundaries
  __attribute__ ((unused)) UInt nb_boundaries= mesh.getNbElementGroups();
  AKANTU_DEBUG_INFO(nb_boundaries << " boundaries advertised initially by Mesh.");  

  AKANTU_DEBUG_INFO("Building boundaries");

  // Two methods: either building using data loaded from the mesh file in MeshData
  // or build with automatic numbering
  mesh.createGroupsFromMeshData<std::string>("physical_names");

  // Second inspection of the number of boundaries (should not be 0)
  nb_boundaries = mesh.getNbElementGroups();

  AKANTU_DEBUG_INFO(nb_boundaries << " boundaries advertised by Mesh.");
  AKANTU_DEBUG_ASSERT(nb_boundaries != 0, "No boundary detected!");

  std::cout << (*dynamic_cast<GroupManager*>(&mesh)) << std::endl;

  akantu::finalize();

  return EXIT_SUCCESS;
}


