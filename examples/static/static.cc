/**
 * @file   static.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Jan 18 2016
 *
 * @brief  This code refers to the implicit static example from the user manual
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
#include "non_linear_solver.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

#define bar_length 0.01
#define bar_height 0.01

/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  initialize("material.dat", argc, argv);

  const UInt spatial_dimension = 2;

  Mesh mesh(spatial_dimension);
  mesh.read("square.msh");

  SolidMechanicsModel model(mesh);

  /// model initialization
  model.initFull(_analysis_method = _static);

  model.setBaseName("static");
  model.addDumpFieldVector("displacement");
  model.addDumpField("external_force");
  model.addDumpField("internal_force");
  model.addDumpField("grad_u");

  /// Dirichlet boundary conditions
  model.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "Fixed_x");
  model.applyBC(BC::Dirichlet::FixedValue(0.0, _y), "Fixed_y");
  model.applyBC(BC::Dirichlet::FixedValue(0.0001, _y), "Traction");
  model.dump();

  auto & solver = model.getNonLinearSolver();
  solver.set("max_iterations", 2);
  solver.set("threshold", 2e-4);
  solver.set("convergence_type", SolveConvergenceCriteria::_solution);

  model.solveStep();

  model.dump();

  finalize();

  return EXIT_SUCCESS;
}
