/**
 * @file   test_heat_transfer_model_cube3d.cc
 *
 * @author Srinivasa Babu Ramisetti <srinivasa.ramisetti@epfl.ch>
 * @author Rui Wang <rui.wang@epfl.ch>
 *
 * @date creation: Sun May 01 2011
 * @date last modification: Fri Jan 26 2018
 *
 * @brief  test of the class HeatTransferModel on the 3d cube
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
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "heat_transfer_model.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"

/* -------------------------------------------------------------------------- */
using namespace akantu;

UInt spatial_dimension = 3;
ElementType type = _tetrahedron_4;

/* -------------------------------------------------------------------------- */

int main(int argc, char * argv[]) {
  initialize("material.dat", argc, argv);

  Mesh mesh(spatial_dimension);
  mesh.read("cube.msh");

  HeatTransferModel model(mesh);
  // initialize everything
  model.initFull();

  // get and set stable time step
  Real time_step = model.getStableTimeStep() * 0.8;
  std::cout << "Stable Time Step is : " << time_step / .8 << std::endl;
  std::cout << "time step is:" << time_step << std::endl;
  model.setTimeStep(time_step);

  /// boundary conditions
  const Array<Real> & nodes = mesh.getNodes();
  Array<bool> & boundary = model.getBlockedDOFs();
  Array<Real> & temperature = model.getTemperature();
  UInt nb_nodes = mesh.getNbNodes();

  double length;
  length = 1.;

  for (UInt i = 0; i < nb_nodes; ++i) {
    temperature(i) = 100.;

    // to insert a heat source
    Real dx = nodes(i, 0) - length / 2.;
    Real dy = nodes(i, 1) - length / 2.;
    Real dz = nodes(i, 2) - length / 2.;
    Real d = sqrt(dx * dx + dy * dy + dz * dz);
    if (d < 0.1) {
      boundary(i) = true;
      temperature(i) = 300.;
    }
  }

  model.setBaseName("heat_transfer_cube3d");
  model.addDumpField("temperature");
  model.addDumpField("temperature_rate");
  model.addDumpField("internal_heat_rate");

  // //for testing
  int max_steps = 1000;

  for (int i = 0; i < max_steps; i++) {
    model.solveStep();

    if (i % 100 == 0)
      model.dump();

    if (i % 10 == 0) {
      std::cout << "Step " << i << "/" << max_steps << std::endl;
    }
  }

  return 0;
}
