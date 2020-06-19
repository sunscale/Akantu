/**
 * @file   test_material_iterative_strength_reduction.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Thu Nov 26 12:20:15 2015
 *
 * @brief  test the material iterative stiffness reduction
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
#include "material_damage_iterative.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */
using namespace akantu;
/* -------------------------------------------------------------------------- */
/* Main                                                                       */
/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  Math::setTolerance(1e-13);
  debug::setDebugLevel(dblWarning);

  initialize("material_stiffness_reduction.dat", argc, argv);

  const UInt spatial_dimension = 2;
  ElementType element_type = _triangle_3;
  const auto & comm = Communicator::getStaticCommunicator();
  Int prank = comm.whoAmI();

  /// read the mesh and partion it
  Mesh mesh(spatial_dimension);
  if (prank == 0) {
    mesh.read("two_elements.msh");
  }

  mesh.distribute();
  /// model creation
  SolidMechanicsModel model(mesh);
  /// initialization of the model
  model.initFull(SolidMechanicsModelOptions(_static));

  /// boundary conditions
  /// Dirichlet BC
  mesh.createGroupsFromMeshData<std::string>(
      "physical_names"); // creates groups from mesh names
  model.applyBC(BC::Dirichlet::FixedValue(0, _x), "left");
  model.applyBC(BC::Dirichlet::FixedValue(0, _y), "bottom");
  model.applyBC(BC::Dirichlet::FixedValue(2., _y), "top");

  /// add fields that should be dumped
  model.setBaseName("material_iterative_stiffness_reduction_test");
  model.addDumpField("material_index");
  model.addDumpFieldVector("displacement");
  ;
  model.addDumpField("stress");
  model.addDumpField("blocked_dofs");
  model.addDumpField("residual");
  model.addDumpField("grad_u");
  model.addDumpField("damage");
  model.addDumpField("partitions");
  model.addDumpField("Sc");
  model.addDumpField("force");
  model.addDumpField("equivalent_stress");
  model.addDumpField("ultimate_strain");

  model.dump();

  MaterialDamageIterative<spatial_dimension> & material =
      dynamic_cast<MaterialDamageIterative<spatial_dimension> &>(
          model.getMaterial(0));

  UInt nb_damaged_elements = 0;
  Real E = material.get("E");
  std::cout << std::setprecision(12);
  const Array<Real> & damage =
      material.getInternal<Real>("damage")(element_type, _not_ghost);
  const Array<Real> & Sc =
      material.getInternal<Real>("Sc")(element_type, _not_ghost);

  /// solve the system
  do {
    model.solveStep();

    nb_damaged_elements = material.updateDamage();

    for (UInt e = 0; e < mesh.getNbElement(element_type, _not_ghost); ++e) {
      std::cout << "the new modulus is " << (1 - damage(0)) * E << std::endl;
      std::cout << "the new strength is " << Sc(0) << std::endl;
    }
    model.dump();

  } while (nb_damaged_elements);

  finalize();

  return EXIT_SUCCESS;
}
