/**
 * @file   explicit_dynamic.cc
 *
 * @author Seyedeh Mohadeseh Taheri Mousavi <mohadeseh.taherimousavi@epfl.ch>
 *
 * @date creation: Sun Jul 12 2015
 * @date last modification: Mon Jan 18 2016
 *
 * @brief  This code refers to the explicit dynamic example from the user manual
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
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */
#include <fstream>
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char * argv[]) {
  initialize("material.dat", argc, argv);

  const UInt spatial_dimension = 3;
  const Real pulse_width = 2.;
  const Real A = 0.01;
  Real time_step;
  Real time_factor = 0.8;
  UInt max_steps = 1000;

  Mesh mesh(spatial_dimension);

  if (Communicator::getStaticCommunicator().whoAmI() == 0)
    mesh.read("bar.msh");

  SolidMechanicsModel model(mesh);

  /// model initialization
  model.initFull(_analysis_method = _explicit_lumped_mass);

  time_step = model.getStableTimeStep();
  std::cout << "Time Step = " << time_step * time_factor << "s (" << time_step
            << "s)" << std::endl;
  model.setTimeStep(time_step * time_factor);

  /// boundary and initial conditions
  Array<Real> & displacement = model.getDisplacement();
  const Array<Real> & nodes = mesh.getNodes();

  for (UInt n = 0; n < mesh.getNbNodes(); ++n) {
    Real x = nodes(n) - 2;
    // Sinus * Gaussian
    Real L = pulse_width;
    Real k = 0.1 * 2 * M_PI * 3 / L;
    displacement(n) = A * sin(k * x) * exp(-(k * x) * (k * x) / (L * L));
  }

  std::ofstream energy;
  energy.open("energy.csv");

  energy << "id,rtime,epot,ekin,tot" << std::endl;

  model.setBaseName("explicit_dynamic");
  model.addDumpField("displacement");
  model.addDumpField("velocity");
  model.addDumpField("acceleration");
  model.addDumpField("stress");
  model.dump();

  for (UInt s = 1; s <= max_steps; ++s) {
    model.solveStep();

    Real epot = model.getEnergy("potential");
    Real ekin = model.getEnergy("kinetic");

    energy << s << "," << s * time_step << "," << epot << "," << ekin << ","
           << epot + ekin << "," << std::endl;

    if (s % 10 == 0)
      std::cout << "passing step " << s << "/" << max_steps << std::endl;
    model.dump();
  }

  energy.close();

  finalize();

  return EXIT_SUCCESS;
}
