/**
 * @file   test_solid_mechanics_model_bar_traction2d_parallel.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
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
#include "mpi.h"
#include <fstream>
#include <limits>

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"
#include "communicator.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char * argv[]) {
  UInt spatial_dimension = 2;
  UInt max_steps = 5000;
  Real time_factor = 0.8;

  initialize("material.dat", argc, argv);

  Mesh mesh(spatial_dimension);

  const auto & comm = Communicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  debug::setDebugLevel(dblWarning);

  if (prank == 0)
    mesh.read("bar2.msh");

  mesh.distribute();

  /* ------------------------------------------------------------------------ */
  /* Initialization                                                           */
  /* ------------------------------------------------------------------------ */
  SolidMechanicsModel model(mesh);
  model.initFull();

  if (prank == 0)
    std::cout << model.getMaterial(0) << std::endl;

  model.assembleMassLumped();

  UInt nb_nodes = mesh.getNbNodes();

  /* ------------------------------------------------------------------------ */
  /* Boundary + initial conditions                                            */
  /* ------------------------------------------------------------------------ */
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

  model.synchronizeBoundaries();

  model.setBaseName("bar2d_parallel");
  model.addDumpField("displacement");
  model.addDumpField("mass");
  model.addDumpField("velocity");
  model.addDumpField("acceleration");
  model.addDumpFieldVector("external_force");
  model.addDumpFieldVector("internal_force");
  model.addDumpField("stress");
  model.addDumpField("strain");
  model.addDumpField("partitions");
  model.dump();

  std::ofstream energy;
  if (prank == 0) {
    energy.open("energy_bar_2d_para.csv");
    energy << "id,rtime,epot,ekin,tot" << std::endl;
  }

  double total_time = 0.;

  /// Setting time step
  Real time_step = model.getStableTimeStep() * time_factor;
  if (prank == 0)
    std::cout << "Time Step = " << time_step << "s" << std::endl;
  model.setTimeStep(time_step);

  /* ------------------------------------------------------------------------ */
  /* Main loop                                                                */
  /* ------------------------------------------------------------------------ */
  for (UInt s = 1; s <= max_steps; ++s) {

    double start = MPI_Wtime();

    model.solveStep();

    double end = MPI_Wtime();

    Real epot = model.getEnergy("potential");
    Real ekin = model.getEnergy("kinetic");

    total_time += end - start;

    if (prank == 0 && (s % 100 == 0)) {
      std::cerr << "passing step " << s << "/" << max_steps << std::endl;
    }
    energy << s << "," << (s - 1) * time_step << "," << epot << "," << ekin
           << "," << epot + ekin << std::endl;

    if (s % 100 == 0) {
      model.dump();
    }
  }

  if (prank == 0)
    std::cout << "Time : " << psize << " " << total_time / max_steps << " "
              << total_time << std::endl;

  finalize();
  return EXIT_SUCCESS;
}
