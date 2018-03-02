/**
 * @file   test_heat_transfer_model_square2d_implicit.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sun May 01 2011
 * @date last modification: Mon Jan 29 2018
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
#include <string>
/* -------------------------------------------------------------------------- */

using namespace akantu;

/* -------------------------------------------------------------------------- */
UInt spatial_dimension = 2;
std::string base_name;

int main(int argc, char * argv[]) {
  initialize("material.dat", argc, argv);

  // create mesh
  Mesh mesh(spatial_dimension);
  mesh.read("square.msh");

  HeatTransferModel model(mesh);
  // initialize everything
  model.initFull(_analysis_method = _static);

  // boundary conditions
  const Array<Real> & nodes = mesh.getNodes();
  Array<bool> & blocked_dofs = model.getBlockedDOFs();
  Array<Real> & temperature = model.getTemperature();
  double length;
  Real dx, dy, dz;
  length = 1.;
  UInt nb_nodes = nodes.size();
  for (UInt i = 0; i < nb_nodes; ++i) {
    temperature(i) = 100.;

    dx = nodes(i, 0) - length / 4.;
    if (spatial_dimension > 1)
      dy = nodes(i, 1) - length / 4.;
    if (spatial_dimension == 3)
      dz = nodes(i, 2) - length / 4.;
    Real d = sqrt(dx * dx + dy * dy + dz * dz);
    //    if(dx < 0.0){
    if (d < 0.1) {
      blocked_dofs(i) = true;
      temperature(i) = 300.;
    }
  }

  // model.assembleInternalHeatRate();
  model.setBaseName("heat_transfer_static_2d");
  model.addDumpField("temperature");
  model.addDumpField("internal_heat_rate");
  model.addDumpField("conductivity");
  model.addDumpField("blocked_dofs");
  model.dump();

  model.solveStep();
  model.dump();

  return 0;
}
