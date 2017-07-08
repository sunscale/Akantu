/**
 * @file   test_weight_computation.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 * @date creation: Sat Sep 26 2015
 * @date last modification: Wed Oct 07 2015
 *
 * @brief  test for the weight computation with base weight function
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
#include "dumper_paraview.hh"
#include "non_local_manager.hh"
#include "non_local_neighborhood.hh"
#include "solid_mechanics_model.hh"
#include "test_material.hh"
/* -------------------------------------------------------------------------- */
using namespace akantu;
/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  akantu::initialize("material_weight_computation.dat", argc, argv);

  // some configuration variables
  const UInt spatial_dimension = 2;

  // mesh creation and read
  Mesh mesh(spatial_dimension);
  mesh.read("plate.msh");

  /// model creation
  SolidMechanicsModel model(mesh);

  /// model initialization changed to use our material
  model.initFull();

  /// dump material index in paraview
  model.addDumpField("material_index");
  model.dump();

  /// save the weights in a file
  NonLocalNeighborhood<BaseWeightFunction> & neighborhood =
    dynamic_cast<NonLocalNeighborhood<BaseWeightFunction> &>(
        model.getNeighborhood("test_region"));

  neighborhood.saveWeights("weights");
  /// print results to screen for validation
  std::ifstream weights;
  weights.open("weights.0");
  std::string current_line;
  while (getline(weights, current_line))
    std::cout << current_line << std::endl;
  weights.close();
  finalize();

  return EXIT_SUCCESS;
}
