/**
 * @file   test_structural_mechanics_model_kirchhoff_shell_patch_test_4_5_5.cc
 *
 * @author Damien Spielmann <damien.spielmann@epfl.ch>
 *
 * @date creation: Thu Feb 21 2013
 * @date last modification: Wed Nov 22 2017
 *
 * @brief  patch test exemple 4.5.5 c.f. modelisation des structures par
 * éléments finis J.-L. Batoz/G Dhatt
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
#include <fstream>
#include <limits>

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "material.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "structural_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

#define TYPE _discrete_kirchhoff_triangle_18

using namespace akantu;

int main(int argc, char * argv[]) {
  initialize(argc, argv);
  Mesh shell(3);
  debug::setDebugLevel(dblWarning);
  std::cout << "Initialisation" << std::endl;
  /* --------------------------------------------------------------------------
   */
  // Defining the mesh

  UInt nb_nodes = 5;
  UInt nb_element = 4;

  Array<Real> & nodes = const_cast<Array<Real> &>(shell.getNodes());
  nodes.resize(nb_nodes);

  Real a = 20.;
  Real b = 10.;

  nodes(0, 0) = 0.;
  nodes(0, 1) = 0.;
  nodes(0, 2) = 0.;

  nodes(1, 0) = 2 * a;
  nodes(1, 1) = 0.;
  nodes(1, 2) = 0.;

  nodes(2, 0) = 0.;
  nodes(2, 1) = 2 * b;
  nodes(2, 2) = 0.;

  nodes(3, 0) = 2 * a;
  nodes(3, 1) = 2 * b;
  nodes(3, 2) = 0.;

  nodes(4, 0) = 15.;
  nodes(4, 1) = 15.;
  nodes(4, 2) = 0.;

  shell.addConnectivityType(TYPE);
  Array<UInt> & connectivity =
      const_cast<Array<UInt> &>(shell.getConnectivity(TYPE));

  connectivity.resize(nb_element);

  connectivity(0, 0) = 1;
  connectivity(0, 1) = 3;
  connectivity(0, 2) = 4;

  connectivity(1, 0) = 3;
  connectivity(1, 1) = 2;
  connectivity(1, 2) = 4;

  connectivity(2, 0) = 2;
  connectivity(2, 1) = 4;
  connectivity(2, 2) = 0;

  connectivity(3, 0) = 0;
  connectivity(3, 1) = 1;
  connectivity(3, 2) = 4;

  akantu::MeshIOMSH mesh_io;
  mesh_io.write("b_beam_3_12_10_13.msh", shell);
  std::cout << "Mesh definition" << std::endl;
  /* --------------------------------------------------------------------------
   */
  // Defining the materials

  akantu::StructuralMechanicsModel model(shell); // ä döfinir

  StructuralMaterial mat1;
  mat1.E = 1000;
  mat1.nu = 0.3;
  mat1.t = 1;

  model.addMaterial(mat1);

  std::cout << "Material Definition" << std::endl;
  /* --------------------------------------------------------------------------
   */
  // Defining the deplacement

  model.initFull();

  // Array<Real> & forces = model.getForce();
  Array<Real> & displacement = model.getDisplacement();
  Array<bool> & boundary = model.getBlockedDOFs();

  displacement(0, 0) = 0;
  displacement(0, 1) = 0;
  displacement(0, 2) = 0;
  displacement(0, 3) = 0;
  displacement(0, 4) = 0;

  displacement(1, 0) = 0;
  displacement(1, 1) = 0;
  displacement(1, 2) = -800;
  displacement(1, 3) = -40;
  displacement(1, 4) = -20;

  displacement(2, 0) = 0;
  displacement(2, 1) = 0;
  displacement(2, 2) = -200;
  displacement(2, 3) = -10;
  displacement(2, 4) = -20;

  displacement(3, 0) = 0;
  displacement(3, 1) = 0;
  displacement(3, 2) = -1400;
  displacement(3, 3) = -50;
  displacement(3, 4) = -40;

  /*displacement(4,0)=0;
    displacement(4,1)=0;*/
  /* displacement(4,2)=;
  displacement(4,3)=;
  displacement(4,4)=;*/

  /* --------------------------------------------------------------------------
   */
  // Defining the boundary conditions

  boundary(0, 0) = true;
  boundary(0, 1) = true;
  boundary(0, 2) = true;
  boundary(0, 3) = true;
  boundary(0, 4) = true;
  boundary(0, 5) = true;

  boundary(1, 0) = true;
  boundary(1, 1) = true;
  boundary(1, 2) = true;
  boundary(1, 3) = true;
  boundary(1, 4) = true;
  boundary(1, 5) = true;

  boundary(2, 0) = true;
  boundary(2, 1) = true;
  boundary(2, 2) = true;
  boundary(2, 3) = true;
  boundary(2, 4) = true;
  boundary(2, 5) = true;

  boundary(3, 0) = true;
  boundary(3, 1) = true;
  boundary(3, 2) = true;
  boundary(3, 3) = true;
  boundary(3, 4) = true;
  boundary(3, 5) = true;

  // boundary(4,0) = true;
  // boundary(4,1) = true;
  // boundary(4,2) = true;
  // boundary(4,3) = true;
  // boundary(4,4) = true;
  boundary(4, 5) = true;

  std::cout << "BC Definition" << std::endl;
  /* --------------------------------------------------------------------------
   */
  // Solve

  Real error;

  model.assembleStiffnessMatrix();
  std::cout << "Assemble Done" << std::endl;
  model.getStiffnessMatrix().saveMatrix("K_4_5_5.mtx");

  UInt count = 0;
  std::cout << "Matrix saved" << std::endl;

  model.addDumpField("displacement");
  model.addDumpField("rotation");
  model.addDumpField("force");
  model.addDumpField("momentum");

  do {
    model.updateResidual();
    model.solve();
    count++;
  } while (!model.testConvergenceIncrement(1e-10, error) && count < 10);

  /* --------------------------------------------------------------------------
   */
  // Post-Processing

  model.computeStresses();

  // const SparseMatrix = model.getStiffnessMatrix();
  std::cout << "u  = " << displacement(4, 0) << std::endl;
  std::cout << "v  = " << displacement(4, 1) << std::endl;
  std::cout << "w5  = " << displacement(4, 2) << std::endl;
  std::cout << "betax  = " << displacement(4, 3) << std::endl;
  std::cout << "betay  = " << displacement(4, 4) << std::endl;
  std::cout << "betaz  = " << displacement(4, 5) << std::endl;

  // model.dump();
}
