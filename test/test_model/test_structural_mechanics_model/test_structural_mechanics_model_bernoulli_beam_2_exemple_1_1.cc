/**
 * @file   test_structural_mechanics_model_bernoulli_beam_2_exemple_1_1.cc
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 *
 * @date creation: Fri Jul 15 2011
 * @date last modification: Sun Oct 19 2014
 *
 * @brief  Computation of the analytical exemple 1.1 in the TGC vol 6
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
#include "aka_common.hh"
#include "mesh.hh"
#include "sparse_matrix_aij.hh"
#include "structural_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char * argv[]) {
  initialize(argc, argv);
  constexpr ElementType type = _bernoulli_beam_2;
  constexpr UInt dim = 2;

  Mesh mesh(dim);

  // Pushing nodes
  auto & nodes = mesh.getNodes();
  Vector<Real> node = {0, 0};
  nodes.push_back(node);
  node = {10, 0};
  nodes.push_back(node);
  node = {18, 0};
  nodes.push_back(node);

  // Pushing connectivity
  mesh.addConnectivityType(type);
  auto & connectivity = mesh.getConnectivity(type);
  Vector<UInt> element = {0, 1};
  connectivity.push_back(element);
  element = {1, 2};
  connectivity.push_back(element);

  StructuralMechanicsModel model(mesh);

  StructuralMaterial mat;
  mat.E = 3e10;
  mat.I = 0.0025;
  mat.A = 0.01;

  model.addMaterial(mat);

  mat.E = 3e10;
  mat.I = 0.00128;
  mat.A = 0.01;

  model.addMaterial(mat);
  // mat.E = 1;
  // mat.I = 1;
  // mat.A = 1;
  // model.addMaterial(mat);
  // model.addMaterial(mat);

  model.initFull();

  // Boundary conditions
  auto boundary = model.getBlockedDOFs().begin(3);
  // clang-format off
  *boundary = {true, true, true}; ++boundary;
  *boundary = {false, true, false}; ++boundary;
  *boundary = {false, true, false}; ++boundary;
  // clang-format on

  // Forces
  Real M = 3600;  // Nm
  Real q = -6000; // kN/m
  Real L = 10;    // m
  auto & forces = model.getExternalForce();
  forces(2, 2) = -M; // moment on last node
  forces(0, 1) = q * L / 2;
  forces(0, 2) = q * L * L / 12;
  forces(1, 1) = q * L / 2;
  forces(1, 2) = -q * L * L / 12;
  forces(2, 0) = mat.E * mat.A / 18;

  // Materials
  auto & materials = model.getElementMaterial(type);
  materials(0) = 0;
  materials(1) = 1;

  model.solveStep();

  // dynamic_cast<const SparseMatrixAIJ &>(model.getDOFManager().getMatrix("K"))
  //     .saveMatrix("stiffness.mtx");

  auto d1 = model.getDisplacement()(1, 2);
  auto d2 = model.getDisplacement()(2, 2);
  auto d3 = model.getDisplacement()(1, 0);

  if (!Math::are_float_equal(d1, 5.6 / 4800) ||  // first rotation
      !Math::are_float_equal(d2, -3.7 / 4800) || // second rotation
      !Math::are_float_equal(d3, 10 / 18.))      // axial deformation
    return 1;

  return 0;
}
