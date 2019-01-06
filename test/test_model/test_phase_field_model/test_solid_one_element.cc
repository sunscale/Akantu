/**
 * @file   tets_phase_field_2d.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Mon Oct 1 2018
 *
 * @brief  test of the class PhaseFieldModel on the 2d square
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
#include "non_linear_solver.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;
const UInt spatial_dimension = 2;
/* -------------------------------------------------------------------------- */
void applyDisplacement(SolidMechanicsModel &);
/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[]) {

  initialize("material_solid.dat", argc, argv);
  
  Mesh mesh(spatial_dimension);
  mesh.read("test_one_element.msh");

  SolidMechanicsModel model(mesh);
  model.initFull(_analysis_method = _static);

  model.setBaseName("solid_one_element");
  model.addDumpField("stress");
  model.addDumpField("grad_u");
  model.addDumpField("displacement");
  model.dump();

  UInt nbSteps = 1000;

  for (UInt s = 0; s < nbSteps; ++s) {
    applyDisplacement(model);
    model.solveStep();
    model.dump();
  }
  
  finalize();

  return EXIT_SUCCESS;

}


/* -------------------------------------------------------------------------- */
void applyDisplacement(SolidMechanicsModel & model) {
  auto & displacement = model.getDisplacement();

  auto & positions = model.getMesh().getNodes();
  auto & blocked_dofs = model.getBlockedDOFs();

  
  for (UInt n = 0; n < model.getMesh().getNbNodes(); ++n) {
    if (positions(n, 1) == -0.5) {
      displacement(n, 0) = 0;
      displacement(n, 1) = 0;
      blocked_dofs(n, 0) = true;
      blocked_dofs(n ,1) = true;
    }
    else {
      displacement(n, 0) = 0;
      displacement(n, 1) += 1.e-4;
      blocked_dofs(n, 0) = true;
      blocked_dofs(n ,1) = true;
    }
  }
}
