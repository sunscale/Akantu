/**
 * @file   test_mesh_boundary.cc
 *
 * @author Dana Christen <dana.christen@gmail.com>
 *
 * @date creation: Fri May 03 2013
 * @date last modification: Sun Aug 13 2017
 *
 * @brief  Thest the element groups
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "mesh.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char * argv[]) {
  UInt spatialDimension(3);

  akantu::initialize(argc, argv);

  Mesh mesh(spatialDimension, "mesh_names");

  std::cout << "Loading the mesh." << std::endl;

  mesh.read("./cube_physical_names.msh");

  std::cout << "Examining mesh:" << std::endl;

  // Inspection of the number of boundaries
  UInt nb_boundaries = mesh.getNbElementGroups(spatialDimension - 1);
  AKANTU_DEBUG_INFO(nb_boundaries << " boundaries advertised by Mesh.");
  if (nb_boundaries == 0) {
    std::cout << "No boundary detected!" << std::endl;
    return 1;
  }

  std::cout << (*dynamic_cast<GroupManager *>(&mesh)) << std::endl;

  akantu::finalize();

  return 0;
}
