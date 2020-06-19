/**
 * @file   test_periodic_plate.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Thu Jan 21 10:11:04 2016
 *
 * @brief  Test for correct application of periodic boundary conditions
 *
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

/* -------------------------------------------------------------------------- */
#include <iostream>

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model_RVE.hh"

using namespace akantu;

int main(int argc, char * argv[]) {

  akantu::initialize("material_test_boundary.dat", argc, argv);

  const UInt spatial_dimension = 2;
  Mesh mesh(spatial_dimension);
  mesh.read("periodic_plate.msh");

  SolidMechanicsModelRVE model(mesh, false);
  auto mat_selector = std::make_shared<MeshDataMaterialSelector<std::string>>(
      "physical_names", model);
  model.setMaterialSelector(mat_selector);

  /// model initialization
  model.initFull();

  /// apply macroscopic deformation gradient at corner nodes
  /// consider a constant strain field
  Matrix<Real> grad_u_macro(spatial_dimension, spatial_dimension, 0.);
  grad_u_macro(0, 1) = 1.;
  model.applyBoundaryConditions(grad_u_macro);

  model.setBaseName("periodic-plate");
  model.addDumpFieldVector("displacement");
  model.addDumpField("stress");
  model.addDumpField("grad_u");
  model.addDumpField("blocked_dofs");
  model.addDumpField("material_index");
  // model.addDumpField      (""      );
  model.dump();

  model.solveStep();

  Real average_strain = model.averageTensorField(0, 1, "strain");
  std::cout << "the average strain is: " << average_strain << std::endl;

  model.dump();

  finalize();

  return EXIT_SUCCESS;
}
