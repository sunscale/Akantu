/**
 * @file   parallel_2d.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Jan 18 2016
 *
 * @brief  Parallel example
 *
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
#include "communicator.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char * argv[]) {
  initialize("material.dat", argc, argv);

  UInt spatial_dimension = 2;
  UInt max_steps = 10000;
  Real time_factor = 0.8;
  Real max_disp = 1e-6;

  Mesh mesh(spatial_dimension);

  const auto & comm = Communicator::getStaticCommunicator();
  Int prank = comm.whoAmI();
  if (prank == 0) {
    // Read the mesh
    mesh.read("square_2d.msh");
  }

  mesh.distribute();

  SolidMechanicsModel model(mesh);
  model.initFull();

  if (prank == 0)
    std::cout << model.getMaterial(0) << std::endl;

  model.setBaseName("multi");
  model.addDumpFieldVector("displacement");
  model.addDumpFieldVector("velocity");
  model.addDumpFieldVector("acceleration");
  model.addDumpFieldTensor("stress");
  model.addDumpFieldTensor("grad_u");

  /// boundary conditions
  Real eps = 1e-16;
  const Array<Real> & pos = mesh.getNodes();
  Array<Real> & disp = model.getDisplacement();
  Array<bool> & boun = model.getBlockedDOFs();

  Real left_side = mesh.getLowerBounds()(0);
  Real right_side = mesh.getUpperBounds()(0);

  for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
    if (std::abs(pos(i, 0) - left_side) < eps) {
      disp(i, 0) = max_disp;
      boun(i, 0) = true;
    }

    if (std::abs(pos(i, 0) - right_side) < eps) {
      disp(i, 0) = -max_disp;
      boun(i, 0) = true;
    }
  }

  Real time_step = model.getStableTimeStep() * time_factor;
  std::cout << "Time Step = " << time_step << "s" << std::endl;
  model.setTimeStep(time_step);

  model.dump();
  for (UInt s = 1; s <= max_steps; ++s) {
    model.solveStep();
    if (s % 200 == 0)
      model.dump();

    if (prank == 0 && s % 100 == 0)
      std::cout << "passing step " << s << "/" << max_steps << std::endl;
  }

  finalize();
  return EXIT_SUCCESS;
}
