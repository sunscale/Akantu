/**
 * @file   embedded.cc
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Tue Dec 01 2015
 * @date last modification: Mon Jan 18 2016
 *
 * @brief  This code gives an example of a simulation using the embedded model
 *
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
#include "embedded_interface_model.hh"
#include "non_linear_solver.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char * argv[]) {
  initialize("material.dat", argc, argv);

  const UInt dim = 2;

  // Loading the concrete mesh
  Mesh mesh(dim);
  mesh.read("concrete.msh");

  // Loading the reinforcement mesh
  Mesh reinforcement_mesh(dim, "reinforcement_mesh");

  // Exception is raised because reinforcement
  // mesh contains only segments, i.e. 1D elements
  try {
    reinforcement_mesh.read("reinforcement.msh");
  } catch (debug::Exception & e) {
  }

  // Model creation
  EmbeddedInterfaceModel model(mesh, reinforcement_mesh, dim);
  model.initFull(EmbeddedInterfaceModelOptions(_static));

  // Boundary conditions
  model.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "XBlocked");
  model.applyBC(BC::Dirichlet::FixedValue(0.0, _y), "YBlocked");

  Vector<Real> force(dim);
  force(0) = 0.0;
  force(1) = -1.0;

  model.applyBC(BC::Neumann::FromTraction(force), "Force");

  // Dumping the concrete
  model.setBaseName("concrete");
  model.addDumpFieldVector("displacement");
  model.addDumpFieldVector("external_force");
  model.addDumpFieldVector("internal_force");
  model.addDumpFieldTensor("stress");

  // Dumping the reinforcement
  model.setBaseNameToDumper("reinforcement", "reinforcement");
  model.addDumpFieldTensorToDumper(
      "reinforcement", "stress_embedded"); // dumping stress in reinforcement

  auto & solver = model.getNonLinearSolver();
  solver.set("max_iterations", 1);
  solver.set("threshold", 1e-6);
  solver.set("convergence_type", SolveConvergenceCriteria::_residual);

  model.solveStep();

  // Dumping model
  model.dump();
  model.dump("reinforcement");

  finalize();
  return EXIT_SUCCESS;
}
