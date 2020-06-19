/**
 * @file   test_material_standard_linear_solid_deviatoric_relaxation_tension.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Aug 09 2010
 * @date last modification: Mon Jun 12 2017
 *
 * @brief  test of the viscoelastic material: relaxation
 *
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
#include "solid_mechanics_model.hh"

using namespace akantu;

int main(int argc, char * argv[]) {
  akantu::initialize("material_standard_linear_solid_deviatoric_relaxation.dat",
                     argc, argv);

  // sim data
  Real T = 10.;
  Real eps = 0.001;

  //  const UInt dim = 3;
  const UInt dim = 2;
  Real sim_time = 25.;
  // Real sim_time = 250.;
  Real time_factor = 0.1;

  Real tolerance = 1e-5;

  Mesh mesh(dim);
  mesh.read("test_material_standard_linear_solid_deviatoric_relaxation.msh");
  //  mesh_io.read("hexa_structured.msh",mesh);
  // const ElementType element_type = _hexahedron_8;
  const ElementType element_type = _quadrangle_4;
  SolidMechanicsModel model(mesh);

  /* ------------------------------------------------------------------------ */
  /* Initialization                                                           */
  /* ------------------------------------------------------------------------ */
  model.initFull();
  std::cout << model.getMaterial(0) << std::endl;

  model.assembleMassLumped();

  model.assembleInternalForces();
  model.getMaterial(0).setToSteadyState();

  std::stringstream filename_sstr;
  filename_sstr << "test_material_standard_linear_solid_deviatoric_relaxation_"
                   "tension.out";
  std::ofstream output_data;
  output_data.open(filename_sstr.str().c_str());
  output_data << "#[1]-time [2]-sigma_analytic [3+]-sigma_measurements"
              << std::endl;

  Material & mat = model.getMaterial(0);
  const Array<Real> & stress = mat.getStress(element_type);

  Real Eta = mat.get("Eta");
  Real EV = mat.get("Ev");
  Real Einf = mat.get("Einf");
  Real E0 = mat.get("E");

  Real kpa = mat.get("kapa");
  Real mu = mat.get("mu");

  Real gamma = EV / E0;
  Real gammainf = Einf / E0;

  Real tau = Eta / EV;
  std::cout << "relaxation time = " << tau << std::endl;

  UInt nb_nodes = mesh.getNbNodes();
  const Array<Real> & coordinate = mesh.getNodes();
  Array<Real> & displacement = model.getDisplacement();

  /// Setting time step
  Real time_step = model.getStableTimeStep() * time_factor;
  std::cout << "Time Step = " << time_step << "s" << std::endl;
  model.setTimeStep(time_step);

  UInt max_steps = sim_time / time_step;
  UInt out_interval = 1;

  Real time = 0.;

  /* ------------------------------------------------------------------------ */
  /* Main loop                                                                */
  /* ------------------------------------------------------------------------ */
  for (UInt s = 0; s <= max_steps; ++s) {

    if (s % 1000 == 0)
      std::cerr << "passing step " << s << "/" << max_steps << std::endl;

    time = s * time_step;
    // impose displacement
    Real epsilon = 0.;
    if (time < T) {
      epsilon = eps * time / T;
    } else {
      epsilon = eps;
    }

    for (UInt n = 0; n < nb_nodes; ++n) {
      for (UInt d = 0; d < dim; ++d)
        displacement(n, d) = epsilon * coordinate(n, d);
    }

    // compute stress
    model.assembleInternalForces();

    // print output
    if (s % out_interval == 0) {
      // analytical solution
      Real epskk = dim * eps;
      Real solution = 0.;
      if (time < T) {
        solution =
            2 * mu * (eps - epskk / 3.) / T *
                (gammainf * time + gamma * tau * (1 - exp(-time / tau))) +
            gammainf * kpa * epskk * time / T;
      } else {
        solution =
            2 * mu * (eps - epskk / 3.) *
                (gammainf +
                 gamma * tau / T * (exp((T - time) / tau) - exp(-time / tau))) +
            gammainf * kpa * epskk;
      }
      output_data << s * time_step << " " << solution;

      // data output
      Array<Real>::const_matrix_iterator stress_it = stress.begin(dim, dim);
      Array<Real>::const_matrix_iterator stress_end = stress.end(dim, dim);
      for (; stress_it != stress_end; ++stress_it) {
        output_data << " " << (*stress_it)(1, 1);

        // test error
        Real rel_error_1 = std::abs(((*stress_it)(1, 1) - solution) / solution);
        if (rel_error_1 > tolerance) {
          std::cerr << "Relative error: " << rel_error_1 << std::endl;
          return EXIT_FAILURE;
        }
      }
      output_data << std::endl;
    }
  }

  finalize();

  std::cout << "Test successful!" << std::endl;
  return EXIT_SUCCESS;
}
