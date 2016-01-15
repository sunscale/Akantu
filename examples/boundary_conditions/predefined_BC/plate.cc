/**
 * @file   plate.cc
 *
 * @author Dana Christen <dana.christen@epfl.ch>
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Dec 16 2015
 *
 * @brief  boundary condition example
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
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char *argv[]) {
  initialize("material.dat", argc, argv);

  //Do I need to define the element type?
  UInt spatial_dimension = 2;
  UInt max_steps = 3000;

  Mesh mesh(spatial_dimension);
  mesh.read("square.msh");
  mesh.createGroupsFromMeshData<std::string>("physical_names");

  SolidMechanicsModel model(mesh);

  /// model initialization
  model.initFull();

  std::cout << model.getMaterial(0) << std::endl;
  model.assembleMassLumped();

  /// Dirichlet boundary conditions
  model.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "Fixed_x");
  model.applyBC(BC::Dirichlet::FixedValue(0.0, _y), "Fixed_y");

  model.updateResidual();

  model.setBaseName("plate");
  model.addDumpFieldVector("displacement");
  model.addDumpField("blocked_dofs");
  model.addDumpField("force"       );
  model.addDumpField("residual"    );
  model.dump();

  Real time_step = model.getStableTimeStep() * .8;
  model.setTimeStep(time_step);

  Vector<Real> traction(2);
  for(UInt s = 1; s <= max_steps; ++s) {
    if(s % 100 == 0)
      std::cout << "passing step " << s << "/" << max_steps << std::endl;

    // Neumann boundary condition
    traction(0) = s*0.1;
    traction(1) = 0.0;
    model.applyBC(BC::Neumann::FromTraction(traction), "Traction");

    model.explicitPred();
    model.updateResidual();
    model.updateAcceleration();
    model.explicitCorr();

    if(s % 10 == 0)
      model.dump();
  }

  finalize();
  return EXIT_SUCCESS;
}
