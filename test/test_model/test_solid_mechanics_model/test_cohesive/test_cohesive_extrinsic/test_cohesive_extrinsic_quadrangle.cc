/**
 * @file   test_cohesive_extrinsic_quadrangle.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Tue May 08 2012
 * @date last modification: Thu Dec 14 2017
 *
 * @brief  Test for extrinsic cohesive elements and quadrangles
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
#include <limits>

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char * argv[]) {
  initialize("material.dat", argc, argv);

  const UInt spatial_dimension = 2;
  const UInt max_steps = 1000;

  Mesh mesh(spatial_dimension);
  mesh.read("quadrangle.msh");

  SolidMechanicsModelCohesive model(mesh);

  /// model initialization
  model.initFull(
      SolidMechanicsModelCohesiveOptions(_explicit_lumped_mass, true));
  model.getElementInserter().setLimit(_y, -0.05, 0.05);
  model.updateAutomaticInsertion();

  Real time_step = model.getStableTimeStep() * 0.05;
  model.setTimeStep(time_step);
  //  std::cout << "Time step: " << time_step << std::endl;

  model.assembleMassLumped();

  const Array<Real> & position = mesh.getNodes();
  Array<Real> & velocity = model.getVelocity();
  Array<bool> & boundary = model.getBlockedDOFs();
  Array<Real> & displacement = model.getDisplacement();
  //  const Array<Real> & residual = model.getResidual();

  UInt nb_nodes = mesh.getNbNodes();

  /// boundary conditions
  for (UInt n = 0; n < nb_nodes; ++n) {
    if (position(n, 1) > 0.99 || position(n, 1) < -0.99)
      boundary(n, 1) = true;

    if (position(n, 0) > 0.99 || position(n, 0) < -0.99)
      boundary(n, 0) = true;
  }

  model.assembleInternalForces();

  /// initial conditions
  Real loading_rate = 0.2;
  Real disp_update = loading_rate * time_step;
  for (UInt n = 0; n < nb_nodes; ++n) {
    velocity(n, 1) = loading_rate * position(n, 1);
  }

  /// Main loop
  for (UInt s = 1; s <= max_steps; ++s) {

    /// update displacement on extreme nodes
    for (UInt n = 0; n < nb_nodes; ++n) {
      if (position(n, 1) > 0.99 || position(n, 1) < -0.99)
        displacement(n, 1) += disp_update * position(n, 1);
    }

    model.checkCohesiveStress();

    model.solveStep();

    if (s % 1 == 0) {
      std::cout << "passing step " << s << "/" << max_steps << std::endl;
    }
  }

  mesh.write("mesh_final.msh");

  Real Ed = model.getEnergy("dissipated");

  Real Edt = 200;

  std::cout << Ed << " " << Edt << std::endl;

  if (Ed < Edt * 0.99 || Ed > Edt * 1.01 || std::isnan(Ed)) {
    std::cout << "The dissipated energy is incorrect" << std::endl;
    return EXIT_FAILURE;
  }

  finalize();

  std::cout << "OK: test_cohesive_extrinsic_quadrangle was passed!"
            << std::endl;
  return EXIT_SUCCESS;
}
