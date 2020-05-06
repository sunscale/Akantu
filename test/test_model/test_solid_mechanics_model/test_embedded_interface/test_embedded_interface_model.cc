/**
 * @file   test_embedded_interface_model.cc
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Wed Mar 25 2015
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  Embedded model test based on potential energy
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

#include <iostream>

#include "aka_common.hh"
#include "embedded_interface_model.hh"
#include "sparse_matrix.hh"

using namespace akantu;

int main(int argc, char * argv[]) {
  debug::setDebugLevel(dblWarning);
  initialize("material.dat", argc, argv);

  UInt dim = 2;
  Math::setTolerance(1e-7);

  // Mesh here is a 1x1 patch
  Mesh mesh(dim);
  mesh.read("embedded_mesh.msh");

  Array<Real> nodes_vec(2, dim, "reinforcement_nodes");
  nodes_vec.storage()[0] = 0;
  nodes_vec.storage()[1] = 0.5;
  nodes_vec.storage()[2] = 1;
  nodes_vec.storage()[3] = 0.5;

  Array<UInt> conn_vec(1, 2, "reinforcement_connectivity");
  conn_vec.storage()[0] = 0;
  conn_vec.storage()[1] = 1;

  Array<std::string> names_vec(1, 1, "reinforcement", "reinforcement_names");

  Mesh reinforcement_mesh(dim, "reinforcement_mesh");
  reinforcement_mesh.getNodes().copy(nodes_vec);
  reinforcement_mesh.addConnectivityType(_segment_2);
  reinforcement_mesh.getConnectivity(_segment_2).copy(conn_vec);
  reinforcement_mesh.getElementalData<std::string>("physical_names")
      .alloc(1, 1, _segment_2);
  reinforcement_mesh.getData<std::string>("physical_names")(_segment_2)
      .copy(names_vec);

  EmbeddedInterfaceModel model(mesh, reinforcement_mesh, dim);
  model.initFull(_analysis_method = _static);

  Array<Real> & nodes = mesh.getNodes();
  Array<Real> & forces = model.getExternalForce();
  Array<bool> & bound = model.getBlockedDOFs();

  forces(2, 0) = -250;
  forces(5, 0) = -500;
  forces(8, 0) = -250;

  for (UInt i = 0; i < mesh.getNbNodes(); i++) {
    if (Math::are_float_equal(nodes(i, 0), 0.))
      bound(i, 0) = true;
    if (Math::are_float_equal(nodes(i, 1), 0.))
      bound(i, 1) = true;
  }

  model.addDumpFieldVector("displacement");
  model.addDumpFieldTensor("stress");

  model.setBaseNameToDumper("reinforcement", "reinforcement");
  model.addDumpFieldTensorToDumper("reinforcement", "stress_embedded");

  model.solveStep();

  model.getDOFManager().getMatrix("K").saveMatrix("matrix_test");

  model.dump();

  Real pot_energy = model.getEnergy("potential");

  if (std::abs(pot_energy - 7.37343e-06) > 1e-5)
    return EXIT_FAILURE;

  finalize();
  return 0;
}
