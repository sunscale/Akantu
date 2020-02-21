/**
 * @file   test_remove_damage_weight_function.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Oct 07 2015
 * @date last modification: Tue Dec 05 2017
 *
 * @brief  Test the damage weight funcion for non local computations
 *
 * @section LICENSE
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
#include "dumper_paraview.hh"
#include "non_local_manager.hh"
#include "non_local_neighborhood.hh"
#include "solid_mechanics_model.hh"
#include "test_material.hh"
/* -------------------------------------------------------------------------- */
using namespace akantu;
/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  akantu::initialize("material_remove_damage.dat", argc, argv);

  // some configuration variables
  const UInt spatial_dimension = 2;
  ElementType element_type = _quadrangle_4;
  GhostType ghost_type = _not_ghost;

  // mesh creation and read
  Mesh mesh(spatial_dimension);
  mesh.read("plate.msh");

  /// model creation
  SolidMechanicsModel model(mesh);
  /// creation of material selector
  auto && mat_selector =
      std::make_shared<MeshDataMaterialSelector<std::string>>("physical_names",
                                                              model);
  model.setMaterialSelector(mat_selector);

  /// model initialization changed to use our material
  model.initFull();
  /// dump material index in paraview
  model.addDumpField("material_index");
  model.addDumpField("grad_u");
  model.addDumpField("grad_u non local");
  model.addDumpField("damage");
  model.dump();

  /// apply constant strain field in all elements except element 3 and 15
  Matrix<Real> applied_strain(spatial_dimension, spatial_dimension);
  applied_strain.clear();
  for (UInt i = 0; i < spatial_dimension; ++i)
    applied_strain(i, i) = 2.;

  /// apply different strain in element 3 and 15
  Matrix<Real> modified_strain(spatial_dimension, spatial_dimension);
  modified_strain.clear();
  for (UInt i = 0; i < spatial_dimension; ++i)
    modified_strain(i, i) = 1.;

  /// apply constant grad_u field in all elements
  for (UInt m = 0; m < model.getNbMaterials(); ++m) {
    Material & mat = model.getMaterial(m);
    Array<Real> & grad_u = const_cast<Array<Real> &>(
        mat.getInternal<Real>("eigen_grad_u")(element_type, ghost_type));

    auto grad_u_it = grad_u.begin(spatial_dimension, spatial_dimension);
    auto grad_u_end = grad_u.end(spatial_dimension, spatial_dimension);
    UInt element_counter = 0;
    for (; grad_u_it != grad_u_end; ++grad_u_it, ++element_counter)
      if (element_counter == 12 || element_counter == 13 ||
          element_counter == 14 || element_counter == 15)
        (*grad_u_it) = -1. * modified_strain;
      else
        (*grad_u_it) = -1. * applied_strain;
  }

  /// compute the non-local strains
  model.assembleInternalForces();
  model.dump();
  /// save the weights in a file
  auto & neighborhood_1 = model.getNonLocalManager().getNeighborhood("mat_1");
  auto & neighborhood_2 = model.getNonLocalManager().getNeighborhood("mat_2");
  neighborhood_1.saveWeights("before_0");
  neighborhood_2.saveWeights("before_1");
  for (UInt n = 0; n < 2; ++n) {
    /// print results to screen for validation
    std::stringstream sstr;
    sstr << "before_" << n << ".0";
    std::ifstream weights;
    weights.open(sstr.str());
    std::string current_line;
    while (getline(weights, current_line))
      std::cout << current_line << std::endl;
    weights.close();
  }

  /// apply damage to not have the elements with lower strain impact the
  /// averaging
  for (UInt m = 0; m < model.getNbMaterials(); ++m) {
    auto & mat =
        dynamic_cast<MaterialDamage<spatial_dimension> &>(model.getMaterial(m));

    auto & damage = const_cast<Array<Real> &>(
        mat.getInternal<Real>("damage")(element_type, ghost_type));

    auto dam_it = damage.begin();
    auto dam_end = damage.end();
    UInt element_counter = 0;
    for (; dam_it != dam_end; ++dam_it, ++element_counter)
      if (element_counter == 12 || element_counter == 13 ||
          element_counter == 14 || element_counter == 15)
        *dam_it = 0.9;
  }

  /// compute the non-local strains
  model.assembleInternalForces();
  neighborhood_1.saveWeights("after_0");
  neighborhood_2.saveWeights("after_1");
  for (UInt n = 0; n < 2; ++n) {
    /// print results to screen for validation
    std::stringstream sstr;
    sstr << "after_" << n << ".0";
    std::ifstream weights;
    weights.open(sstr.str());
    std::string current_line;
    while (getline(weights, current_line))
      std::cout << current_line << std::endl;
    weights.close();
  }

  model.dump();

  /// verify the result: non-local averaging over constant field must
  /// yield same constant field
  Real test_result = 0.;
  Matrix<Real> difference(spatial_dimension, spatial_dimension, 0.);
  Matrix<Real> difference_in_damaged_elements(spatial_dimension,
                                              spatial_dimension, 0.);
  for (UInt m = 0; m < model.getNbMaterials(); ++m) {
    difference_in_damaged_elements.clear();
    auto & mat = model.getMaterial(m);
    auto & grad_u_nl = const_cast<Array<Real> &>(
        mat.getInternal<Real>("grad_u non local")(element_type, ghost_type));

    auto grad_u_nl_it = grad_u_nl.begin(spatial_dimension, spatial_dimension);
    auto grad_u_nl_end = grad_u_nl.end(spatial_dimension, spatial_dimension);
    UInt element_counter = 0;
    for (; grad_u_nl_it != grad_u_nl_end; ++grad_u_nl_it, ++element_counter) {
      if (element_counter == 12 || element_counter == 13 ||
          element_counter == 14 || element_counter == 15)
        difference_in_damaged_elements += (*grad_u_nl_it);
      else
        difference = (*grad_u_nl_it) - applied_strain;
      test_result += difference.norm<L_2>();
    }
    difference_in_damaged_elements *= (1 / 4.);
    difference_in_damaged_elements -= (1.41142 * modified_strain);
    test_result += difference_in_damaged_elements.norm<L_2>();
  }

  if (test_result > 10.e-5) {
    std::cout << "the total norm is: " << test_result << std::endl;
    return EXIT_FAILURE;
  }

  finalize();

  return EXIT_SUCCESS;
}
