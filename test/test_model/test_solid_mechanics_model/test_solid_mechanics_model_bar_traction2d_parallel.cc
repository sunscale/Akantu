/**
 * @file   test_solid_mechanics_model_bar_traction2d_parallel.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sun Sep 12 2010
 * @date last modification: Thu Jun 05 2014
 *
 * @brief  test of the class SolidMechanicsModel
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
#include <limits>
#include <fstream>
#include "mpi.h"

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"

using namespace akantu;

int main(int argc, char *argv[]) {
  akantu::UInt spatial_dimension = 2;
  akantu::UInt max_steps = 5000;
  akantu::Real time_factor = 0.8;

  akantu::initialize("material.dat", argc, argv);

  akantu::Mesh mesh(spatial_dimension);

  akantu::StaticCommunicator & comm =
    akantu::StaticCommunicator::getStaticCommunicator();
  akantu::Int psize = comm.getNbProc();
  akantu::Int prank = comm.whoAmI();

  akantu::debug::setDebugLevel(akantu::dblWarning);

  akantu::MeshPartition * partition = NULL;
  if(prank == 0) {
    mesh.read("bar2.msh");
    partition = new akantu::MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);
  }

  /* ------------------------------------------------------------------------ */
  /* Initialization                                                           */
  /* ------------------------------------------------------------------------ */
  akantu::SolidMechanicsModel model(mesh);

  /// model initialization
  model.initParallel(partition);

  model.initFull();

  if(prank == 0)
    std::cout << model.getMaterial(0) << std::endl;

  model.assembleMassLumped();


  akantu::UInt nb_nodes = mesh.getNbNodes();

  /* ------------------------------------------------------------------------ */
  /* Boundary + initial conditions                                            */
  /* ------------------------------------------------------------------------ */
  akantu::Real eps = 1e-16;
  const akantu::Array<akantu::Real> & pos = mesh.getNodes();
  akantu::Array<akantu::Real> & disp = model.getDisplacement();
  akantu::Array<bool> & boun = model.getBlockedDOFs();

  for (akantu::UInt i = 0; i < nb_nodes; ++i) {
    if(pos(i, 0) >= 9.) disp(i, 0) = (pos(i, 0) - 9) / 100.;
    if(pos(i) <= eps)   boun(i, 0) = true;
    if(pos(i, 1) <= eps || pos(i, 1) >= 1 - eps ) boun(i, 1) = true;
  }

  model.synchronizeBoundaries();

  model.updateResidual();

  model.setBaseName("bar2d_parallel");
  model.addDumpField("displacement");
  model.addDumpField("mass"        );
  model.addDumpField("velocity"    );
  model.addDumpField("acceleration");
  model.addDumpField("force"       );
  model.addDumpField("residual"    );
  model.addDumpField("stress"      );
  model.addDumpField("strain"      );
  model.addDumpField("partitions"  );
  model.dump();

  std::ofstream energy;
  if(prank == 0) {
    energy.open("energy_bar_2d_para.csv");
    energy << "id,rtime,epot,ekin,tot" << std::endl;
  }

  double total_time = 0.;

  /// Setting time step
  akantu::Real time_step = model.getStableTimeStep() * time_factor;
  if(prank == 0)
    std::cout << "Time Step = " << time_step << "s" << std::endl;
  model.setTimeStep(time_step);

  /* ------------------------------------------------------------------------ */
  /* Main loop                                                                */
  /* ------------------------------------------------------------------------ */
  for(akantu::UInt s = 1; s <= max_steps; ++s) {

    double start = MPI_Wtime();

    model.explicitPred();

    model.updateResidual();
    model.updateAcceleration();
    model.explicitCorr();

    double end = MPI_Wtime();

    akantu::Real epot = model.getPotentialEnergy();
    akantu::Real ekin = model.getKineticEnergy();

    total_time += end - start;

    if(prank == 0 && (s % 100 == 0)) {
      std::cerr << "passing step " << s << "/" << max_steps << std::endl;
    }
    energy << s << "," << (s-1)*time_step << "," << epot << "," << ekin << "," << epot + ekin
	     << std::endl;

    if(s % 100 == 0) {
      model.dump();
    }
  }

  if(prank == 0) std::cout << "Time : " << psize << " " << total_time / max_steps << " " << total_time << std::endl;

  akantu::finalize();
  return EXIT_SUCCESS;
}
