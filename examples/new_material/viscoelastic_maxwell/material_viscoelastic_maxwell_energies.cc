/**
 * @file   material_viscoelastic_maxwell_energies.cc
 *
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 *
 * @date creation: Tue Nov 20 2018
 * @date last modification:
 *
 * @brief  Example of using viscoelastic material and computing energies
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
/* -------------------------------------------------------------------------- */
#include "material_viscoelastic_maxwell.hh"
#include "non_linear_solver.hh"
#include "solid_mechanics_model.hh"
#include "sparse_matrix.hh"

using namespace akantu;

/* -------------------------------------------------------------------------- */
/* Main                                                                       */
/* -------------------------------------------------------------------------- */

int main(int argc, char * argv[]) {
  akantu::initialize("material_viscoelastic_maxwell.dat", argc, argv);

  // sim data
  Real eps = 0.1;

  const UInt dim = 2;
  Real sim_time = 100.;
  Real T = 10.;
  Mesh mesh(dim);
  mesh.read("material_viscoelastic_maxwell_mesh.msh");

  SolidMechanicsModel model(mesh);

  /* ------------------------------------------------------------------------ */
  /* Initialization                                                           */
  /* ------------------------------------------------------------------------ */
  model.initFull(_analysis_method = _static);
  std::cout << model.getMaterial(0) << std::endl;

  std::stringstream filename_sstr;
  filename_sstr << "material_viscoelastic_maxwell_output.out";
  std::ofstream output_data;
  output_data.open(filename_sstr.str().c_str());

  Material & mat = model.getMaterial(0);

  Real time_step = 0.1;

  UInt nb_nodes = mesh.getNbNodes();
  const Array<Real> & coordinate = mesh.getNodes();
  Array<Real> & displacement = model.getDisplacement();
  Array<bool> & blocked = model.getBlockedDOFs();

  /// Setting time step

  model.setTimeStep(time_step);

  model.setBaseName("dynamic");
  model.addDumpFieldVector("displacement");
  model.addDumpField("blocked_dofs");
  model.addDumpField("external_force");
  model.addDumpField("internal_force");
  model.addDumpField("grad_u");
  model.addDumpField("stress");
  model.addDumpField("strain");

  UInt max_steps = sim_time / time_step + 1;
  Real time = 0.;

  auto & solver = model.getNonLinearSolver();
  solver.set("max_iterations", 10);
  solver.set("threshold", 1e-7);
  solver.set("convergence_type", SolveConvergenceCriteria::_residual);

  /* ------------------------------------------------------------------------ */
  /* Main loop                                                                */
  /* ------------------------------------------------------------------------ */
  for (UInt s = 0; s <= max_steps; ++s) {

    std::cout << "Time Step = " << time_step << "s" << std::endl;
    std::cout << "Time = " << time << std::endl;

    // impose displacement
    Real epsilon = 0;
    if (time < T) {
      epsilon = eps * time / T;
    } else {
      epsilon = eps;
    }

    for (UInt n = 0; n < nb_nodes; ++n) {
      if (Math::are_float_equal(coordinate(n, 0), 0.0)) {
        displacement(n, 0) = 0;
        blocked(n, 0) = true;
        displacement(n, 1) = epsilon * coordinate(n, 1);
        blocked(n, 1) = true;
      } else if (Math::are_float_equal(coordinate(n, 1), 0.0)) {
        displacement(n, 0) = epsilon * coordinate(n, 0);
        blocked(n, 0) = true;
        displacement(n, 1) = 0;
        blocked(n, 1) = true;
      } else if (Math::are_float_equal(coordinate(n, 0), 0.001)) {
        displacement(n, 0) = epsilon * coordinate(n, 0);
        blocked(n, 0) = true;
        displacement(n, 1) = epsilon * coordinate(n, 1);
        blocked(n, 1) = true;
      } else if (Math::are_float_equal(coordinate(n, 1), 0.001)) {
        displacement(n, 0) = epsilon * coordinate(n, 0);
        blocked(n, 0) = true;
        displacement(n, 1) = epsilon * coordinate(n, 1);
        blocked(n, 1) = true;
      }
    }

    try {
      model.solveStep();
    } catch (debug::Exception & e) {
    }

    // for debugging
    // auto int_force = model.getInternalForce();
    // auto &K = model.getDOFManager().getMatrix("K");
    // K.saveMatrix("K.mtx");

    Int nb_iter = solver.get("nb_iterations");
    Real error = solver.get("error");
    bool converged = solver.get("converged");

    if (converged) {
      std::cout << "Converged in " << nb_iter << " iterations" << std::endl;
    } else {
      std::cout << "Didn't converge after " << nb_iter
                << " iterations. Error is " << error << std::endl;
      return EXIT_FAILURE;
    }

    model.dump();

    Real epot = mat.getEnergy("potential");
    Real edis = mat.getEnergy("dissipated");
    Real work = mat.getEnergy("work");

    // data output
    output_data << s * time_step << " " << epsilon << " " << epot << " " << edis
                << " " << work << std::endl;
    time += time_step;
  }
  output_data.close();
  finalize();
}
