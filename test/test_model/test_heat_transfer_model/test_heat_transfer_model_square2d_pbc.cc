/**
 * @file   test_heat_transfer_model_square2d_pbc.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Sun May 01 2011
 * @date last modification: Wed Feb 03 2016
 *
 * @brief  test of the class HeatTransferModel on the 3d cube
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
#include "heat_transfer_model.hh"
/* -------------------------------------------------------------------------- */
#include <fstream>
#include <iostream>
#include <string.h>
using namespace std;
/* -------------------------------------------------------------------------- */

int main(int argc, char * argv[]) {
  akantu::UInt spatial_dimension = 2;
  akantu::initialize("material.dat", argc, argv);

  // create mesh
  akantu::Mesh mesh(spatial_dimension);
  mesh.read("square_tri3.msh");

  akantu::HeatTransferModel model(mesh);
  // initialize everything
  model.initFull();

  // initialize PBC
  model.setPBC(1, 1, 1);
  model.initPBC();

  // assemble the lumped capacity
  model.assembleCapacityLumped();

  // get stable time step
  akantu::Real time_step = model.getStableTimeStep() * 0.8;
  cout << "time step is:" << time_step << endl;
  model.setTimeStep(time_step);

  // boundary conditions
  const akantu::Array<akantu::Real> & nodes =
      model.getFEEngine().getMesh().getNodes();
  akantu::Array<bool> & boundary = model.getBlockedDOFs();
  akantu::Array<akantu::Real> & temperature = model.getTemperature();
  double length;
  length = 1.;
  akantu::UInt nb_nodes = model.getFEEngine().getMesh().getNbNodes();
  for (akantu::UInt i = 0; i < nb_nodes; ++i) {
    temperature(i) = 100.;

    akantu::Real dx = nodes(i, 0) - length / 4.;
    akantu::Real dy = 0.0;
    akantu::Real dz = 0.0;

    if (spatial_dimension > 1)
      dy = nodes(i, 1) - length / 4.;
    if (spatial_dimension == 3)
      dz = nodes(i, 2) - length / 4.;
    akantu::Real d = sqrt(dx * dx + dy * dy + dz * dz);
    //    if(dx < 0.0){
    if (d < 0.1) {
      boundary(i) = true;
      temperature(i) = 300.;
    }
  }

  model.updateResidual();
  model.setBaseName("heat_transfer_square_2d_pbc");
  model.addDumpField("temperature");
  model.addDumpField("temperature_rate");
  model.addDumpField("residual");
  model.addDumpField("capacity_lumped");
  model.dump();

  // main loop
  int max_steps = 1000;
  for (int i = 0; i < max_steps; i++) {
    model.explicitPred();
    model.updateResidual();
    model.solveExplicitLumped();
    model.explicitCorr();
    if (i % 100 == 0)
      model.dump();
    if (i % 10 == 0)
      std::cout << "Step " << i << "/" << max_steps << std::endl;
  }
  cout << "\n\n Stable Time Step is : " << time_step << "\n \n" << endl;

  return 0;
}
