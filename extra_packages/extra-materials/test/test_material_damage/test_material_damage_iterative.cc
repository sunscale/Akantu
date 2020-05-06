/**
 * @file   test_material_damage_iterative.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Thu Nov 26 12:20:15 2015
 *
 * @brief  test the material damage iterative
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

  initialize("material.dat", argc, argv);

  const UInt spatial_dimension = 2;
  ElementType element_type = _triangle_3;
  const auto & comm = Communicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  /// read the mesh and partion it
  Mesh mesh(spatial_dimension);

  if (prank == 0) {
    mesh.read("plate.msh");
  }

  mesh.distribute();

  /// model creation
  SolidMechanicsModel model(mesh);

  /// initialization of the model
  model.initFull(_analysis_method = _static);

  /// boundary conditions
  /// Dirichlet BC
  mesh.createGroupsFromMeshData<std::string>(
      "physical_names"); // creates groups from mesh names
  model.applyBC(BC::Dirichlet::FixedValue(0, _x), "left");
  model.applyBC(BC::Dirichlet::FixedValue(0, _y), "bottom");
  model.applyBC(BC::Dirichlet::FixedValue(2., _y), "top");

  /// add fields that should be dumped
  model.setBaseName("material_damage_iterative_test");
  model.addDumpFieldVector("displacement");
  ;
  model.addDumpField("stress");
  model.addDumpField("blocked_dofs");
  model.addDumpField("residual");
  model.addDumpField("grad_u");
  model.addDumpField("damage");
  model.addDumpField("partitions");
  model.addDumpField("material_index");
  model.addDumpField("Sc");
  model.addDumpField("force");
  model.addDumpField("equivalent_stress");

  model.dump();

  MaterialDamageIterative<spatial_dimension> & material =
      dynamic_cast<MaterialDamageIterative<spatial_dimension> &>(
          model.getMaterial(0));

  UInt nb_damaged_elements = 0;
  Real max_eq_stress = 0;

  /// solve the system
  model.solveStep();

  model.dump();

  /// check that the normalized equivalent stress
  Array<Real> & eq_stress =
      material.getInternal<Real>("equivalent_stress")(element_type, _not_ghost);
  Array<Real>::const_scalar_iterator eq_stress_it = eq_stress.begin();
  UInt nb_elements = mesh.getNbElement(element_type, _not_ghost);
  for (UInt e = 0; e < nb_elements; ++e, ++eq_stress_it) {
    if (!Math::are_float_equal(*eq_stress_it, 0.1)) {
      std::cout << "Error in the equivalent normalized stress" << std::endl;
      finalize();
      return EXIT_FAILURE;
    }
  }

  /// get the maximum equivalent stress
  max_eq_stress = material.getNormMaxEquivalentStress();

  nb_damaged_elements = 0;
  if (max_eq_stress > 1.)
    nb_damaged_elements = material.updateDamage();

  if (nb_damaged_elements) {
    std::cout << "Damage occured even though the normalized stress is below 1"
              << std::endl;
    finalize();
    return EXIT_FAILURE;
  }

  /// weaken material locally to cause damage
  Array<Real> & strength = const_cast<Array<Real> &>(
      material.getInternal<Real>("Sc")(element_type, _not_ghost));
  Array<Real>::scalar_iterator strength_it = strength.begin();
  ++strength_it;
  *strength_it = 0.9;
  strength_it += 4;
  *strength_it = 0.898;

  /// solve the system again
  model.solveStep();

  /// get the maximum equivalent stress
  max_eq_stress = material.getNormMaxEquivalentStress();

  nb_damaged_elements = 0;
  if (max_eq_stress > 1.)
    nb_damaged_elements = material.updateDamage();

  UInt nb_damaged_elements_per_proc = 2;
  if (nb_damaged_elements != psize * nb_damaged_elements_per_proc) {
    std::cout << "Error in number of damaged elements" << std::endl;
    finalize();
    return EXIT_FAILURE;
  }

  /// check that damage occured in correct elements
  Real dam_diff = 0.;
  Array<Real> & damage =
      material.getInternal<Real>("damage")(element_type, _not_ghost);
  Array<Real>::const_scalar_iterator damage_it = damage.begin();
  for (UInt e = 0; e < nb_elements; ++e, ++damage_it) {
    if (e == 1 || e == 5)
      dam_diff += std::abs(0.1 - *damage_it);
    else
      dam_diff += (*damage_it);
  }

  if (dam_diff > 1.e-13) {
    std::cout << "Error in damage pattern" << std::endl;
    finalize();
    return EXIT_FAILURE;
  }

  /// solve to compute the stresses correctly for dumping

  model.solveStep();

  model.dump();

  finalize();

  return EXIT_SUCCESS;
}
