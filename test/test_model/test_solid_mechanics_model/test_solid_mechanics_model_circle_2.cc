/**
 * @file   test_solid_mechanics_model_circle_2.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Aug 04 2010
 * @date last modification: Tue Aug 04 2015
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
#include <iostream>

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */
using namespace akantu;

int main(int argc, char * argv[]) {
  akantu::initialize("material.dat", argc, argv);
  UInt max_steps = 10000;
  Real epot, ekin, wext = 0;

  Mesh mesh(2);
  mesh.read("circle2.msh");
  mesh.createBoundaryGroupFromGeometry();

  SolidMechanicsModel model(mesh);

  /// model initialization
  model.initFull();

  Real time_step = model.getStableTimeStep() / 10.;
  model.setTimeStep(time_step);

  std::cout << "-- Time step : " << time_step << " --" << std::endl;

  model.assembleMassLumped();

  // Boundary condition (Neumann)
  Matrix<Real> stress(2, 2);
  stress.eye(Real(1e3));
  model.applyBC(BC::Neumann::FromHigherDim(stress), "boundary_0");

  model.setBaseName("circle2");
  model.addDumpFieldVector("displacement");
  model.addDumpFieldVector("external_force");
  model.addDumpFieldVector("internal_force");
  model.addDumpField("mass");
  model.addDumpField("velocity");
  model.addDumpField("acceleration");
  model.addDumpField("stress");
  model.addDumpField("strain");
  model.dump();

  for (UInt s = 0; s < max_steps; ++s) {
    model.solveStep();

    epot = model.getEnergy("potential");
    ekin = model.getEnergy("kinetic");
    wext += model.getEnergy("external work");

    std::cout << s << "," << epot << "," << ekin << "," << wext << ","
              << epot + ekin - wext << std::endl;

    if (s % 100 == 0)
      model.dump();
  }

  return EXIT_SUCCESS;
}
