/**
 * @file   test_heat_transfer_model_cube3d_pbc.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Srinivasa Babu Ramisetti <srinivasa.ramisetti@epfl.ch>
 * @author Rui Wang <rui.wang@epfl.ch>
 *
 * @date creation: Sun May 01 2011
 * @date last modification: Thu Jun 05 2014
 *
 * @brief  test of the class HeatTransferModel on the 3d cube
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_common.hh"
#include "mesh.hh"
#include "heat_transfer_model.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
#include <fstream>

using namespace std;
/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[])
{
  akantu::initialize("material.dat", argc, argv);

  akantu::UInt spatial_dimension = 3;
  //create mesh
  akantu::Mesh mesh(spatial_dimension);
  mesh.read("cube_tet4.msh");

  akantu::HeatTransferModel model(mesh);
  //initialize everything
  model.initFull();

  //initialize PBC
  model.setPBC(1,1,1);
  model.initPBC();

  //assemble the lumped capacity
  model.assembleCapacityLumped();

  //get stable time step
  akantu::Real time_step = model.getStableTimeStep()*0.8;
  cout<<"time step is:"<<time_step<<endl;
  model.setTimeStep(time_step);

  // boundary conditions
  const akantu::Array<akantu::Real> & nodes = model.getFEEngine().getMesh().getNodes();
  akantu::Array<bool> & boundary = model.getBlockedDOFs();
  akantu::Array<akantu::Real> & temperature = model.getTemperature();

  //  double t1, t2;
  double length;
  //  t1 = 300.;
  //  t2 = 100.;
  length = 1.;
  akantu::UInt nb_nodes = model.getFEEngine().getMesh().getNbNodes();
  for (akantu::UInt i = 0; i < nb_nodes; ++i) {
    temperature(i) = 100.;
    akantu::Real dx = nodes(i,0) - length/4.;
    akantu::Real dy = nodes(i,1) - length/4.;
    akantu::Real dz = nodes(i,2) - length/4.;

    akantu::Real d = sqrt(dx*dx + dy*dy + dz*dz);
    if(d < 0.1){
      boundary(i) = true;
      temperature(i) = 300.;
    }
  }

  model.updateResidual();
  model.setBaseName("heat_transfer_cube3d_pbc");
  model.addDumpField("temperature"     );
  model.addDumpField("temperature_rate");
  model.addDumpField("residual"        );
  model.addDumpField("capacity_lumped" );
  model.dump();

  /* ------------------------------------------------------------------------ */
  // //for testing
  int max_steps = 1000;
  /* ------------------------------------------------------------------------ */
  for(int i=0; i<max_steps; i++)
    {
      model.explicitPred();
      model.updateResidual();
      model.solveExplicitLumped();
      model.explicitCorr();

      if(i % 100 == 0) model.dump();

      if(i % 10 == 0)
        std::cout << "Step " << i << "/" << max_steps << std::endl;
    }
  cout<< "\n\n Stable Time Step is : " << time_step << "\n \n" <<endl;

  return 0;
}
