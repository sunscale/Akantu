/**
 * @file   test_structural_mechanics_model_bernoulli_beam_3.cc
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Mon Jan 22 2018
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
#include "sparse_solver.hh"
#include "structural_mechanics_model.hh"
// #include "test_gtest_utils.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
/* -------------------------------------------------------------------------- */

using namespace akantu;

TEST(TestBernoulliBeam3, TestDisplacements) {
  constexpr ElementType type = _bernoulli_beam_3;
  constexpr UInt dim = 3;
  const UInt ndof = ElementClass<type>::getNbDegreeOfFreedom();
  const Real a = std::sqrt(2) / 2; // cos(pi/4)

  Mesh mesh(dim, "test_bernoulli_beam_3");

  // Pushing nodes
  auto & nodes = mesh.getNodes();
  Vector<Real> node = {0, 0, 0};
  nodes.push_back(node);
  node = {a, a, 0};
  // node = {1, 0, 0};
  nodes.push_back(node);
  node = {a, -a, 0};
  // node = {0, 1, 0};
  nodes.push_back(node);

  // Pushing connectivity
  mesh.addConnectivityType(type);
  auto & connectivity = mesh.getConnectivity(type);
  Vector<UInt> element = {0, 1};
  connectivity.push_back(element);
  element = {0, 2};
  connectivity.push_back(element);

  // Pushing normals
  auto & normals =
      mesh.registerData<Real>("extra_normal").alloc(0, dim, type, _not_ghost);
  Vector<Real> normal = {0, 0, 1};
  normals.push_back(normal);
  normal = {0, 0, 1};
  normals.push_back(normal);

  // Creating model
  StructuralMechanicsModel model(mesh, dim, "test_bernoulli_beam_3");

  // Unit material
  StructuralMaterial mat;
  mat.E = 1;
  mat.Iz = 1;
  mat.Iy = 1;
  mat.A = 1;
  mat.GJ = 1;
  model.addMaterial(mat);

  model.initFull();

  // Boundary conditions (blocking all DOFs of nodes 2 & 3)
  auto boundary = ++model.getBlockedDOFs().begin(ndof);
  // clang-format off
  *boundary = {true, true, true, true, true, true}; ++boundary;
  *boundary = {true, true, true, true, true, true}; ++boundary;
  // clang-format on

  // Forces
  Real P = 1; // N
  auto & forces = model.getExternalForce();
  forces.clear();
  forces(0, 2) = -P; // vertical force on first node

  // Setting same material for all elements
  model.getElementMaterial(type).set(0);

  try {
    model.solveStep();
  } catch (debug::SingularMatrixException & e) {
    std::cerr << e.what() << std::endl;
    e.matrix.saveMatrix("jacobian.mtx");
  }

  model.getDOFManager().getMatrix("K").saveMatrix("stiffness.mtx");
  model.getDOFManager().getMatrix("J").saveMatrix("jacobian.mtx");

  auto vz = model.getDisplacement()(0, 2);
  auto thy = model.getDisplacement()(0, 4);
  auto thx = model.getDisplacement()(0, 3);
  auto thz = model.getDisplacement()(0, 5);

  Real tol = Math::getTolerance();

  EXPECT_NEAR(vz, -5. / 48, tol);
  EXPECT_NEAR(thy, -std::sqrt(2) / 8, tol);
  EXPECT_NEAR(thz, 0, tol);
  EXPECT_NEAR(thx, 0, tol);
}
