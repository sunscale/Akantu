/**
 * @file   test_material_non_local.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Aug 31 2011
 * @date last modification: Thu Oct 15 2015
 *
 * @brief  test of the main part of the non local materials
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "custom_non_local_test_material.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;
/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[]) {
  akantu::initialize("material.dat", argc, argv);

  // some configuration variables
  const UInt spatial_dimension = 2;

  Mesh mesh(spatial_dimension);

  StaticCommunicator & comm = akantu::StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  // mesh creation and read
  MeshPartition * partition;
  if(prank == 0) {
    mesh.read("mesh.msh");

    /// partition the mesh
    partition = new MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);
  } else {
    partition = NULL;
  }

  /// model creation
  SolidMechanicsModel  model(mesh);
  model.initParallel(partition);
  delete partition;

  /// model initialization changed to use our material
  model.initFull(_no_init_materials = true);

  model.registerNewCustomMaterials< CustomNonLocalTestMaterial<spatial_dimension> >("custom_non_local_test_material");
  model.initMaterials();

  CustomNonLocalTestMaterial<spatial_dimension> & mat = dynamic_cast<CustomNonLocalTestMaterial<spatial_dimension> &>(model.getMaterial("test"));

  if(prank == 0) std::cout << mat << std::endl;

  // Setting up the dumpers + first dump
  model.setBaseName("non_local_material");
  model.addDumpFieldVector("displacement");
  model.addDumpFieldVector("force"       );
  model.addDumpField("partitions"   );
  model.addDumpField("stress"      );
  model.addDumpField("stress"      );
  model.addDumpField("local_damage");
  model.addDumpField("damage"      );

  model.updateResidual();
  model.dump();


  //Array<Real> & damage = mat.getArray("local_damage", _quadrangle_4);
  Array<Real> & damage = mat.getArray<Real>("local_damage", _triangle_3);

  RandGenerator<UInt> gen;

  for (UInt i = 0; i < 1; ++i) {
    UInt g = (gen() / Real(RandGenerator<UInt>::max() - RandGenerator<UInt>::min()))  * damage.getSize();
    std::cout << prank << " -> " << g << std::endl;
    damage(g) = 1.;
  }

  model.updateResidual();
  model.dump();

  akantu::finalize();
  return EXIT_SUCCESS;
}
