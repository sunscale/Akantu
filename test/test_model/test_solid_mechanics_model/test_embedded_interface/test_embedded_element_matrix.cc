/**
 * @file   test_embedded_element_matrix.cc
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Wed Mar 25 2015
 * @date last modification: Fri Feb 09 2018
 *
 * @brief  test of the class EmbeddedInterfaceModel
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#include "embedded_interface_model.hh"
#include "sparse_matrix_aij.hh"
#include "sparse_solver.hh"

using namespace akantu;

int main(int argc, char * argv[]) {
  debug::setDebugLevel(dblWarning);
  initialize("embedded_element.dat", argc, argv);

  constexpr UInt dim = 2;
  constexpr ElementType type = _segment_2;
  const Real height = 0.4;

  Mesh mesh(dim);
  mesh.read("triangle.msh");

  Mesh reinforcement_mesh(dim, "reinforcement_mesh");
  auto & nodes = reinforcement_mesh.getNodes();
  nodes.push_back(Vector<Real>({0, height}));
  nodes.push_back(Vector<Real>({1, height}));

  reinforcement_mesh.addConnectivityType(type);
  auto & connectivity = reinforcement_mesh.getConnectivity(type);
  connectivity.push_back(Vector<UInt>({0, 1}));

  Array<std::string> names_vec(1, 1, "reinforcement", "reinforcement_names");
  reinforcement_mesh.getElementalData<std::string>("physical_names")
      .alloc(1, 1, type);
  reinforcement_mesh.getData<std::string>("physical_names")(type).copy(
      names_vec);

  EmbeddedInterfaceModel model(mesh, reinforcement_mesh, dim);

  model.initFull(_analysis_method = _static);

  if (model.getInterfaceMesh().getNbElement(type) != 1)
    return EXIT_FAILURE;

  if (model.getInterfaceMesh().getSpatialDimension() != 2)
    return EXIT_FAILURE;

  try { // matrix should be singular
    model.solveStep();
  } catch (debug::SingularMatrixException & e) {
    std::cerr << "Matrix is singular, relax, everything is fine :)"
              << std::endl;
  } catch (debug::Exception & e) {
    std::cerr << "Unexpceted error: " << e.what() << std::endl;
    throw e;
  }

  SparseMatrixAIJ & K =
      dynamic_cast<SparseMatrixAIJ &>(model.getDOFManager().getMatrix("K"));
  K.saveMatrix("stiffness.mtx");

  Math::setTolerance(1e-8);

  // Testing the assembled stiffness matrix
  if (!Math::are_float_equal(K(0, 0), 1. - height) ||
      !Math::are_float_equal(K(0, 2), height - 1.) ||
      !Math::are_float_equal(K(2, 0), height - 1.) ||
      !Math::are_float_equal(K(2, 2), 1. - height))
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
