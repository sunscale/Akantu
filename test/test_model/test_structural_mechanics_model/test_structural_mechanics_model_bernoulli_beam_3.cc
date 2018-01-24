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
#include "structural_mechanics_model.hh"
#include "sparse_solver.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char * argv[]) {
  initialize(argc, argv);
  constexpr ElementType type = _bernoulli_beam_3;
  constexpr UInt dim = 3;
  const UInt ndof = ElementClass<type>::getNbDegreeOfFreedom();
  const Real a = std::sqrt(2) / 2; // cos(pi/4)

  Mesh mesh(dim);

  // Pushing nodes
  auto & nodes = mesh.getNodes();
  Vector<Real> node = {0, 0, 0};
  nodes.push_back(node);
  node = {1, 0, 0};
  nodes.push_back(node);
  node = {0, 1, 0};
  nodes.push_back(node);

  // Pushing connectivity
  mesh.addConnectivityType(type);
  auto & connectivity = mesh.getConnectivity(type);
  Vector<UInt> element = {0, 1};
  connectivity.push_back(element);
  element = {0, 2};
  connectivity.push_back(element);

  // Pushing normals
  auto & normals = mesh.registerData<Real>("extra_normal")
                       .alloc(0, dim, type, _not_ghost);
  Vector<Real> normal = {0, 0, 1};
  normals.push_back(normal);
  normal = {0, 0, 1};
  normals.push_back(normal);

  // Creating model
  StructuralMechanicsModel model(mesh);

  // Unit material
  StructuralMaterial mat;
  mat.E = 1;
  mat.Iz = 1;
  mat.Iy = 1;
  mat.A = 1;
  mat.GJ = 1;
  model.addMaterial(mat);
  // mat.E = 0.5;
  // mat.Iz = 1;
  // mat.Iy = 1;
  // mat.A = 1;
  // mat.GJ = 1;
  // model.addMaterial(mat);

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
  forces(0, 2) = -P; // vertical force on first node

  // Setting same material for all elements
  model.getElementMaterial(type).set(0);

  try {
    model.solveStep();
  } catch (debug::SingularMatrixException & e) {
    std::cerr << e.what() << std::endl;
    e.matrix.saveMatrix("jacobian.mtx");
    return 1;
  }

  model.getDOFManager().getMatrix("J").saveMatrix("jacobian.mtx");

  auto vz = model.getDisplacement()(0, 2);
  auto th = model.getDisplacement()(0, 4);


  if (!Math::are_float_equal(vz, -5. / 48.) ||        // vertical deflection
      !Math::are_float_equal(th, -std::sqrt(2) / 8.)) { // y rotation
    std::cout << "vz = " << vz << "\n"
	      << "th = " << th << std::endl;
    for (auto val : make_view(model.getDisplacement(), 1))
      std::cout << val << std::endl;
    return 1;
  }

  return 0;
}
