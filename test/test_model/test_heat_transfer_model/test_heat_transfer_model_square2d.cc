/**
 * @file   test_heat_transfer_model_square2d.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Srinivasa Babu Ramisetti <srinivasa.ramisetti@epfl.ch>
 *
 * @date creation: Sun May 01 2011
 * @date last modification: Sun Oct 19 2014
 *
 * @brief  test of the class HeatTransferModel on the 3d cube
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
#include "heat_transfer_model.hh"
/* -------------------------------------------------------------------------- */
#include <fstream>
#include <iostream>
#include <string>

using namespace akantu;

/* -------------------------------------------------------------------------- */
UInt spatial_dimension = 2;

std::string base_name;

int main(int argc, char * argv[]) {
  initialize("material.dat", argc, argv);

  // create mesh
  Mesh mesh(spatial_dimension);
  mesh.read("square_tri3.msh");

  HeatTransferModel model(mesh);
  // initialize everything
  model.initFull();

  // get stable time step
  Real time_step = model.getStableTimeStep() * 0.8;
  std::cout << "Stable Time Step is: " << time_step << std::endl;
  std::cout << "Time step is: " << time_step << std::endl;
  model.setTimeStep(time_step);

  // boundary conditions
  const Array<Real> & nodes = model.getFEEngine().getMesh().getNodes();
  Array<bool> & boundary = model.getBlockedDOFs();
  Array<Real> & temperature = model.getTemperature();
  double length;
  length = 1.;
  UInt nb_nodes = model.getFEEngine().getMesh().getNbNodes();
  for (UInt i = 0; i < nb_nodes; ++i) {
    temperature(i) = 100.;

    Real dx = nodes(i, 0) - length / 4.;
    Real dy = 0.0;
    Real dz = 0.0;

    if (spatial_dimension > 1)
      dy = nodes(i, 1) - length / 4.;
    if (spatial_dimension == 3)
      dz = nodes(i, 2) - length / 4.;
    Real d = sqrt(dx * dx + dy * dy + dz * dz);
    //    if(dx < 0.0){
    if (d < 0.1) {
      boundary(i) = true;
      temperature(i) = 300.;
    }
  }

  model.assembleInternalHeatRate();
  model.setBaseName("heat_transfer_square2d");
  model.addDumpField("temperature");
  model.addDumpField("temperature_rate");
  model.addDumpField("internal_heat_rate");
  model.addDumpField("capacity_lumped");
  model.dump();

  // main loop
  int max_steps = 1500;
  for (int i = 0; i < max_steps; i++) {
    model.solveStep();

    if (i % 100 == 0)
      model.dump();
    if (i % 10 == 0)
      std::cout << "Step " << i << "/" << max_steps << std::endl;
  }

  return 0;
}
