/**
 * @file   test_non_local_averaging.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sat Sep 26 2015
 * @date last modification: Tue Dec 05 2017
 *
 * @brief  test for non-local averaging of strain
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
  akantu::initialize("material_avg.dat", argc, argv);

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
  model.dump();

  /// apply constant strain field everywhere in the plate
  Matrix<Real> applied_strain(spatial_dimension, spatial_dimension);
  applied_strain.clear();
  for (UInt i = 0; i < spatial_dimension; ++i)
    applied_strain(i, i) = 2.;

  /// apply constant grad_u field in all elements
  for (auto & mat : model.getMaterials()) {
    auto & grad_us =
        mat.getInternal<Real>("eigen_grad_u")(element_type, ghost_type);
    for (auto & grad_u :
         make_view(grad_us, spatial_dimension, spatial_dimension)) {
      grad_u = -1. * applied_strain;
    }
  }

  /// compute the non-local strains
  model.assembleInternalForces();
  model.dump();

  /// verify the result: non-local averaging over constant field must
  /// yield same constant field
  Real test_result = 0.;
  Matrix<Real> difference(spatial_dimension, spatial_dimension, 0.);

  for (auto & mat : model.getMaterials()) {
    auto & grad_us_nl =
        mat.getInternal<Real>("grad_u non local")(element_type, ghost_type);
    for (auto & grad_u_nl :
         make_view(grad_us_nl, spatial_dimension, spatial_dimension)) {
      difference = grad_u_nl - applied_strain;
      test_result += difference.norm<L_2>();
    }
  }

  if (test_result > 10.e-13) {
    AKANTU_EXCEPTION("the total norm is: " << test_result);
  }

  return 0;
}
