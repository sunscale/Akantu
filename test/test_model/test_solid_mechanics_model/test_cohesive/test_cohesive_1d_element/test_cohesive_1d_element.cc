/**
 * @file   test_cohesive_1d_element.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Jun 14 2013
 * @date last modification: Mon Jun 23 2014
 *
 * @brief  Test for 1D cohesive elements
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
#include "solid_mechanics_model_cohesive.hh"

/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char *argv[]) {
  initialize("material.dat", argc, argv);

  const UInt max_steps = 2000;
  const Real strain_rate = 4;

  UInt spatial_dimension = 1;
  Mesh mesh(spatial_dimension, "mesh");
  mesh.read("bar.msh");

  SolidMechanicsModelCohesive model(mesh);
  model.initFull(SolidMechanicsModelCohesiveOptions(_explicit_lumped_mass, true));

  Real time_step = model.getStableTimeStep()*0.01;
  model.setTimeStep(time_step);
  std::cout << "Time step: " << time_step << std::endl;

  model.assembleMassLumped();

  mesh.computeBoundingBox();
  Real posx_max = mesh.getUpperBounds()(0);
  Real posx_min = mesh.getLowerBounds()(0);

  /// initial conditions
  Array<Real> & velocity = model.getVelocity();
  const Array<Real> & position = mesh.getNodes();
  UInt nb_nodes = mesh.getNbNodes();

  for (UInt n = 0; n < nb_nodes; ++n)
    velocity(n) = strain_rate * (position(n) - (posx_max + posx_min) /2.);

  /// boundary conditions
  Array<bool> & boundary = model.getBlockedDOFs();
  Array<Real> & displacement = model.getDisplacement();
  Real disp_increment = strain_rate * (posx_max - posx_min) / 2. * time_step;

  for(UInt node = 0; node < mesh.getNbNodes(); ++node) {
    if(Math::are_float_equal(position(node), posx_min)) {      // left side
      boundary(node) = true;
    }
    if(Math::are_float_equal(position(node), posx_max)) {      // right side
      boundary(node) = true;
    }
  }


  model.synchronizeBoundaries();
  model.updateResidual();

  // model.setBaseName("extrinsic_parallel");
  // model.addDumpFieldVector("displacement");
  // model.addDumpField("velocity"    );
  // model.addDumpField("acceleration");
  // model.addDumpField("residual"    );
  // model.addDumpField("stress");
  // model.addDumpField("strain");
  // model.dump();

  for (UInt s = 1; s <= max_steps; ++s) {

    model.checkCohesiveStress();
    model.solveStep();

    UInt nb_cohesive_elements = mesh.getNbElement(_cohesive_1d_2);

    if (s % 10 == 0) {
      std::cout << "passing step " << s << "/" << max_steps
		<< ", number of cohesive elemets:"
		<< nb_cohesive_elements << std::endl;

      //      model.dump();
    }

    /// update external work and boundary conditions
    for (UInt n = 0; n < mesh.getNbNodes(); ++n) {
      if(Math::are_float_equal(position(n), posx_min))      // left side
    	displacement(n) -= disp_increment;

      if(Math::are_float_equal(position(n), posx_max))      // right side
    	displacement(n) += disp_increment;
    }
  }

  Real Ed = model.getEnergy("dissipated");
  Real Edt = 100 * 3;

  std::cout << Ed << std::endl;
  if (Ed < Edt * 0.999 || Ed > Edt * 1.001 || std::isnan(Ed)) {
    std::cout << "The dissipated energy is incorrect" << std::endl;
    finalize();
    return EXIT_FAILURE;
  }

  finalize();
  return EXIT_SUCCESS;
}
