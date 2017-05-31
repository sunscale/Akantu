/**
 * @file   test_solid_mechanics_model_bar_traction2d_structured.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Dec 13 2010
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
#include <limits>
#include <fstream>
#include <iostream>

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char *argv[]) {
  initialize("material.dat", argc, argv);
  UInt spatial_dimension = 2;
  UInt max_steps = 10000;
  Real time_factor = 0.2;

  Real epot, ekin;

  Mesh mesh(spatial_dimension);
  mesh.read("bar_structured1.msh");

  SolidMechanicsModel model(mesh);

  /// model initialization
  model.initFull();

  std::cout << model.getMaterial(0) << std::endl;

  /// boundary conditions
  Real eps = 1e-16;
  for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
    if(mesh.getNodes()(i) >= 9)
      model.getDisplacement()(i) = (model.getFEEngine().getMesh().getNodes()(i) - 9) / 100. ;

    if(mesh.getNodes()(i) <= eps)
	model.getBlockedDOFs()(i) = true;

    if(mesh.getNodes()(i, 1) <= eps ||
       mesh.getNodes()(i, 1) >= 1 - eps ) {
      model.getBlockedDOFs()(i, 1) = true;
    }
  }

  Real time_step = model.getStableTimeStep() * time_factor;
  std::cout << "Time Step = " << time_step << "s" << std::endl;
  model.setTimeStep(time_step);

  std::ofstream energy;
  energy.open("energy_bar_2d_structured.csv");
  energy << "id,epot,ekin,tot" << std::endl;


  for(UInt s = 1; s <= max_steps; ++s) {
    model.solveStep();

    epot = model.getEnergy("potential");
    ekin = model.getEnergy("kinetic");

    std::cerr << "passing step " << s << "/" << max_steps << std::endl;
    energy << s << "," << epot << "," << ekin << "," << epot + ekin
	   << std::endl;
  }

  energy.close();

  return EXIT_SUCCESS;
}
