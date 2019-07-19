/**
 * @file   test_material_viscoelastic_maxwell_relaxation.cc
 *
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 *
 * @date creation: Thu May 17 2018
 * @date last modification:
 *
 * @brief  test of the viscoelastic material: relaxation
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
  Real tolerance = 1e-6;
  Mesh mesh(dim);
  mesh.read("test_material_viscoelastic_maxwell.msh");

  const ElementType element_type = _quadrangle_4;
  SolidMechanicsModel model(mesh);

  /* ------------------------------------------------------------------------ */
  /* Initialization                                                           */
  /* ------------------------------------------------------------------------ */
  model.initFull(_analysis_method = _static);
  std::cout << model.getMaterial(0) << std::endl;

  std::stringstream filename_sstr;
  filename_sstr << "test_material_viscoelastic_maxwell.out";
  std::ofstream output_data;
  output_data.open(filename_sstr.str().c_str());

  Material & mat = model.getMaterial(0);

  const Array<Real> & stress = mat.getStress(element_type);

  Vector<Real> Eta = mat.get("Eta");
  Vector<Real> Ev = mat.get("Ev");
  Real Einf = mat.get("Einf");
  Real nu = mat.get("nu");
  Real lambda = Eta(0) / Ev(0);
  Real pre_mult = 1 / (1 + nu) / (1 - 2 * nu);
  Real time_step = 0.1;

  UInt nb_nodes = mesh.getNbNodes();
  const Array<Real> & coordinate = mesh.getNodes();
  Array<Real> & displacement = model.getDisplacement();
  Array<bool> & blocked = model.getBlockedDOFs();

  /// Setting time step

  model.setTimeStep(time_step);

  UInt max_steps = sim_time / time_step + 1;
  Real time = 0.;

  auto & solver = model.getNonLinearSolver();
  solver.set("max_iterations", 200);
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

    // analytical solution
    Real solution_11 = 0.;
    if (time < T) {
      solution_11 = pre_mult * eps / T *
                    (Einf * time + lambda * Ev(0) * (1 - exp(-time / lambda)));
    } else {
      solution_11 =
          pre_mult * eps *
          (Einf + lambda * Ev(0) / T *
                      (exp((T - time) / lambda) - exp(-time / lambda)));
    }

    // data output
    output_data << s * time_step << " " << epsilon << " " << solution_11;
    Array<Real>::const_matrix_iterator stress_it = stress.begin(dim, dim);
    Array<Real>::const_matrix_iterator stress_end = stress.end(dim, dim);
    for (; stress_it != stress_end; ++stress_it) {
      output_data << " " << (*stress_it)(0, 0);
      // test error
      Real rel_error_11 =
          std::abs(((*stress_it)(0, 0) - solution_11) / solution_11);
      if (rel_error_11 > tolerance) {
        std::cerr << "Relative error: " << rel_error_11 << std::endl;
        return EXIT_FAILURE;
      }
    }
    output_data << std::endl;
    time += time_step;
  }
  output_data.close();
  finalize();

  std::cout << "Test successful!" << std::endl;
  return EXIT_SUCCESS;
}
