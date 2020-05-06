/**
 * @file   test_material_FE2.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Sun Jan 31 12:27:02 2016
 *
 * @brief  test the material FE2
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
#include "communicator.hh"
#include "material_FE2.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

/* -------------------------------------------------------------------------- */
/* Main                                                                       */
/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  debug::setDebugLevel(dblWarning);

  initialize("material.dat", argc, argv);

  const auto & comm = Communicator::getStaticCommunicator();
  Int prank = comm.whoAmI();

  /// input parameters for the simulation
  const UInt spatial_dimension = 2;
  const ParserSection & parser = getUserParser();
  std::string mesh_file = parser.getParameter("mesh_file");
  Matrix<Real> prestrain_increment = parser.getParameter("prestrain_increment");
  UInt total_steps = parser.getParameter("total_steps");

  Mesh mesh(spatial_dimension);

  if (prank == 0) {
    mesh.read(mesh_file);
  }

  mesh.distribute();

  /// model creation
  SolidMechanicsModel model(mesh);

  /// set the material selector
  auto mat_selector = std::make_shared<MaterialSelector>();
  mat_selector->setFallback(3);
  model.setMaterialSelector(mat_selector);

  model.initFull(SolidMechanicsModelOptions(_static));

  /* --------------------------------------------------------------------------
   */
  /// boundary conditions
  mesh.createGroupsFromMeshData<std::string>(
      "physical_names"); // creates groups from mesh names
  model.applyBC(BC::Dirichlet::FixedValue(0, _x), "bottom");
  model.applyBC(BC::Dirichlet::FixedValue(0, _y), "bottom");
  //  model.applyBC(BC::Dirichlet::FixedValue(1.e-2, _y), "top");

  model.setBaseName("macro_mesh");
  model.addDumpFieldVector("displacement");
  model.addDumpField("stress");
  model.addDumpField("grad_u");
  model.addDumpField("eigen_grad_u");
  model.addDumpField("blocked_dofs");
  model.addDumpField("material_index");
  model.addDumpField("material_stiffness");

  model.dump();

  /// solve system
  model.assembleStiffnessMatrix();

  std::cout << "first solve step" << std::endl;
  model.solveStep();

  std::cout << "second solve step" << std::endl;
  model.solveStep();

  std::cout << "finished solve steps" << std::endl;
  /// simulate the advancement of the reaction
  MaterialFE2<spatial_dimension> & mat =
      dynamic_cast<MaterialFE2<spatial_dimension> &>(
          model.getMaterial("FE2_mat"));
  Matrix<Real> current_prestrain(spatial_dimension, spatial_dimension, 0.);
  for (UInt i = 0; i < total_steps; ++i) {
    model.dump();
    current_prestrain += prestrain_increment;
    mat.advanceASR(current_prestrain);
    model.dump();
    /// solve for new displacement at the macro-scale
    model.solveStep();
  }

  model.dump();

  finalize();

  return EXIT_SUCCESS;
}
