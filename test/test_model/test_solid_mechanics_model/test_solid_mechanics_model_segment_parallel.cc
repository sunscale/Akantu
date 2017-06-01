/**
 * @file   test_solid_mechanics_model_segment_parallel.cc
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
#include <limits>
/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"
#include "static_communicator.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char * argv[]) {
  UInt spatial_dimension = 1;
  UInt max_steps = 10000;
  Real time_factor = 0.2;

  initialize("material.dat", argc, argv);

  Mesh mesh(spatial_dimension);
  StaticCommunicator & comm =
      StaticCommunicator::getStaticCommunicator();
    Int prank = comm.whoAmI();

  debug::setDebugLevel(dblInfo);

  if (prank == 0) mesh.read("segment.msh");

  mesh.distribute();

  SolidMechanicsModel model(mesh);

  /// model initialization
  model.initFull();

  std::cout << model.getMaterial(0) << std::endl;

  model.setBaseName("segment_parallel");
  model.addDumpField("displacement");
  model.addDumpField("mass");
  model.addDumpField("velocity");
  model.addDumpField("acceleration");
  model.addDumpField("external_force");
  model.addDumpField("internal_force");
  model.addDumpField("stress");
  model.addDumpField("strain");

  /// boundary conditions
  for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
    model.getDisplacement().storage()[spatial_dimension * i] =
        model.getFEEngine().getMesh().getNodes().storage()[i] / 100.;

    if (model.getFEEngine().getMesh().getNodes().storage()[spatial_dimension *
                                                           i] <= 1e-15)
      model.getBlockedDOFs().storage()[i] = true;
  }

  Real time_step = model.getStableTimeStep() * time_factor;
  std::cout << "Time Step = " << time_step << "s" << std::endl;
  model.setTimeStep(time_step);

  for (UInt s = 1; s <= max_steps; ++s) {
    model.solveStep();

    Real epot = model.getEnergy("potential");
    Real ekin = model.getEnergy("kinetic");

    if (prank == 0) {
      std::cout << s << " " << epot << " " << ekin << " " << epot + ekin
                << std::endl;
    }

    model.dump();
    if (s % 10 == 0)
      std::cerr << "passing step " << s << "/" << max_steps << std::endl;
  }

  finalize();

  return EXIT_SUCCESS;
}
