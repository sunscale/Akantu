/**
 * @file   test_solid_mechanics_model_bar_traction2d.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Aug 09 2010
 * @date last modification: Thu Aug 06 2015
 *
 * @brief  test of the class SolidMechanicsModel
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
#include <fstream>
#include <iostream>
#include <limits>

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char * argv[]) {
  initialize("material.dat", argc, argv);
  UInt max_steps = 5000;
  Real time_factor = 0.8;

  UInt spatial_dimension = 2;

  Mesh mesh(spatial_dimension);
  mesh.read("bar2.msh");

  SolidMechanicsModel model(mesh);

  UInt nb_nodes = mesh.getNbNodes();

  /// model initialization
  model.initFull();

  std::cout << model.getMaterial(0) << std::endl;

  /// boundary conditions
  Real eps = 1e-16;
  const Array<Real> & pos = mesh.getNodes();
  Array<Real> & disp = model.getDisplacement();
  Array<bool> & boun = model.getBlockedDOFs();

  for (UInt i = 0; i < nb_nodes; ++i) {
    if (pos(i, 0) >= 9.)
      disp(i, 0) = (pos(i, 0) - 9) / 100.;
    if (pos(i) <= eps)
      boun(i, 0) = true;
    if (pos(i, 1) <= eps || pos(i, 1) >= 1 - eps)
      boun(i, 1) = true;
  }

  /// set the time step
  Real time_step = model.getStableTimeStep() * time_factor;
  std::cout << "Time Step = " << time_step << "s" << std::endl;
  model.setTimeStep(time_step);

  /// initialize the paraview output
  mesh.setBaseName("bar_traction_2d");
  model.addDumpField("displacement");
  model.addDumpField("mass");
  model.addDumpField("velocity");
  model.addDumpField("acceleration");
  model.addDumpField("force");
  model.addDumpField("residual");
  model.addDumpFieldTensor("stress");
  model.addDumpField("grad_u");

  model.addDumpGroupField("displacement", "Top");
  model.dumpGroup("Top");
  model.dump();

  std::ofstream energy;
  energy.open("energy_bar_2d.csv");
  energy << "id,rtime,epot,ekin,tot" << std::endl;

  for (UInt s = 1; s <= max_steps; ++s) {
    model.solveStep();

    Real epot = model.getEnergy("potential");
    Real ekin = model.getEnergy("kinetic");

    energy << s << "," << (s - 1) * time_step << "," << epot << "," << ekin
           << "," << epot + ekin << std::endl;


#ifdef AKANTU_USE_IOHELPER
    if (s % 100 == 0) {
      model.dump();
      model.dumpGroup();
    }
#endif // AKANTU_USE_IOHELPER
    if (s % 100 == 0)
      std::cout << "passing step " << s << "/" << max_steps << std::endl;
  }

  energy.close();

  finalize();

  return EXIT_SUCCESS;
}
