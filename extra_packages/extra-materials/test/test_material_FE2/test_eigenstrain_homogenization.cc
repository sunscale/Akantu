/**
 * @file   test_eigenstrain_homogenization.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Sun Jan 31 12:27:02 2016
 *
 * @brief  test the eigenstrain homogenization
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
#include "solid_mechanics_model_RVE.hh"

using namespace akantu;

/* -------------------------------------------------------------------------- */
/* Main                                                                       */
/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {

  akantu::initialize("mesoscale_materials.dat", argc, argv);

  const UInt spatial_dimension = 2;
  Mesh mesh(spatial_dimension);
  mesh.read("one_inclusion.msh");
  SolidMechanicsModelRVE model(mesh, false);
  auto mat_selector = std::make_shared<MeshDataMaterialSelector<std::string>>(
      "physical_names", model);
  model.setMaterialSelector(mat_selector);

  /// model initialization
  model.initFull();

  /// apply boundary conditions
  Matrix<Real> grad_u_macro(spatial_dimension, spatial_dimension, 0.);
  model.applyBoundaryConditions(grad_u_macro);

  model.setBaseName("one-inclusion");
  model.addDumpFieldVector("displacement");
  model.addDumpField("stress");
  model.addDumpField("grad_u");
  model.addDumpField("blocked_dofs");
  model.addDumpField("material_index");
  model.addDumpField("eigen_grad_u");
  model.dump();

  /// apply eigenstrain
  Matrix<Real> prestrain(spatial_dimension, spatial_dimension, 0.);
  for (UInt i = 0; i < spatial_dimension; ++i)
    prestrain(i, i) = 0.02;
  model.advanceASR(prestrain);

  model.dump();

  Matrix<Real> macro_strain(spatial_dimension, spatial_dimension, 0.);
  model.homogenizeEigenGradU(macro_strain);

  std::cout << "the average eigen_gradu is " << macro_strain << std::endl;

  Matrix<Real> exact_eigenstrain(spatial_dimension, spatial_dimension, 0.);
  for (UInt i = 0; i < spatial_dimension; ++i)
    exact_eigenstrain(i, i) = 0.00125;

  macro_strain -= exact_eigenstrain;

  if (macro_strain.norm<L_2>() > 1.e-10) {
    std::cout << "the test failed!!" << std::endl;
    finalize();
    return EXIT_FAILURE;
  }

  finalize();
  return EXIT_SUCCESS;
}
