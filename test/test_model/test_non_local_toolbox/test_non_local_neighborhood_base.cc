/**
 * @file   test_non_local_neighborhood_base.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Wed Sep 23 16:30:12 2015
 *
 * @brief  test for the class NonLocalNeighborhoodBase
 *
 * @section LICENSE
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
#include "solid_mechanics_model.hh"
#include "test_material.hh"
#include "non_local_neighborhood_base.hh"
#include "dumper_paraview.hh"
/* -------------------------------------------------------------------------- */
using namespace akantu;
/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[]) {
  akantu::initialize("material.dat", argc, argv);

  // some configuration variables
  const UInt spatial_dimension = 2;

  // mesh creation and read
  Mesh mesh(spatial_dimension);
  mesh.read("plate.msh");

  /// model creation
  SolidMechanicsModel  model(mesh);

  /// creation of material selector
  MeshDataMaterialSelector<std::string> * mat_selector;
  mat_selector = new MeshDataMaterialSelector<std::string>("physical_names", model);
  model.setMaterialSelector(*mat_selector);
 
  /// model initialization changed to use our material
  model.initFull(SolidMechanicsModelOptions(_static, true));
  model.registerNewCustomMaterials< TestMaterial<spatial_dimension> >("test_material");
  model.initMaterials();
  /// dump material index in paraview
  model.addDumpField("material_index");
  model.dump();
  /// non-local neighborhood creation
  Real non_local_radius = 0.5;
  NonLocalNeighborhoodBase neighborhood(model, non_local_radius);

  /// insert the quadrature points in the neighborhood
  for (UInt m = 0; m < model.getNbMaterials(); ++m) {
    TestMaterial<spatial_dimension> & mat = dynamic_cast<TestMaterial<spatial_dimension> & > (model.getMaterial(m));
    mat.insertQuads(neighborhood);
  }

  /// create the pairs of quadrature points
  neighborhood.updatePairList();

  /// save the pair of quadrature points and the coords of all neighbors
  std::string output_1 = "quadrature_pairs";
  std::string output_2 = "neighborhoods";
  neighborhood.savePairs(output_1);
  neighborhood.saveNeighborCoords(output_2);

  /// print results to screen for validation
  std::ifstream quad_pairs;
  quad_pairs.open("quadrature_pairs.0");
  std::string current_line;
  while(getline(quad_pairs, current_line))
    std::cout << current_line << std::endl;
  quad_pairs.close();
  std::ifstream neighborhoods;
  neighborhoods.open("neighborhoods.0");
  while(getline(neighborhoods, current_line))
    std::cout << current_line << std::endl;
  neighborhoods.close();

  finalize();
  
  return EXIT_SUCCESS;
}