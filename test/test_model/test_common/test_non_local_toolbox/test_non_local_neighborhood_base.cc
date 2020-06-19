/**
 * @file   test_non_local_neighborhood_base.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sat Sep 26 2015
 * @date last modification: Mon Sep 11 2017
 *
 * @brief  test for the class NonLocalNeighborhoodBase
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
#include "my_model.hh"
/* -------------------------------------------------------------------------- */
#include <fstream>
/* -------------------------------------------------------------------------- */

using namespace akantu;

/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  akantu::initialize("material.dat", argc, argv);

  // some configuration variables
  const UInt spatial_dimension = 2;

  // mesh creation and read
  Mesh mesh(spatial_dimension);
  mesh.read("plate.msh");

  /// model creation
  MyModel model(mesh, spatial_dimension);
  const auto & manager = model.getNonLocalManager();
  const auto & neighborhood = manager.getNeighborhood("test_region");

  /// save the pair of quadrature points and the coords of all neighbors
  std::string output_1 = "quadrature_pairs";
  std::string output_2 = "neighborhoods";
  neighborhood.savePairs(output_1);
  neighborhood.saveNeighborCoords(output_2);

  /// print results to screen for validation
  std::ifstream quad_pairs;
  quad_pairs.open("quadrature_pairs.0");

  std::string current_line;
  while (getline(quad_pairs, current_line))
    std::cout << current_line << std::endl;

  quad_pairs.close();

  std::ifstream neighborhoods;
  neighborhoods.open("neighborhoods.0");
  while (getline(neighborhoods, current_line))
    std::cout << current_line << std::endl;
  neighborhoods.close();

  finalize();

  return EXIT_SUCCESS;
}
