/**
 * @file   test_cohesive_1d_element.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Jun 14 2013
 * @date last modification: Wed Jan 10 2018
 *
 * @brief  Test for 1D cohesive elements
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char * argv[]) {
  initialize("material.dat", argc, argv);

  const UInt max_steps = 2000;
  const Real strain_rate = 5;

  UInt spatial_dimension = 1;
  Mesh mesh(spatial_dimension, "mesh");
  mesh.read("bar.msh");

  Math::setTolerance(1e-7);
  SolidMechanicsModelCohesive model(mesh);
  model.initFull(_analysis_method = _explicit_lumped_mass,
                 _is_extrinsic = true);

  auto time_step = model.getStableTimeStep() * 0.01;
  model.setTimeStep(time_step);
  std::cout << "Time step: " << time_step << std::endl;

  auto posx_max = mesh.getUpperBounds()(_x);
  auto posx_min = mesh.getLowerBounds()(_x);

  /// initial conditions
  const auto & position = mesh.getNodes();
  auto & velocity = model.getVelocity();
  auto nb_nodes = mesh.getNbNodes();

  for (UInt n = 0; n < nb_nodes; ++n)
    velocity(n) = strain_rate * (position(n) - (posx_max + posx_min) / 2.);

  /// boundary conditions
  model.applyBC(BC::Dirichlet::FlagOnly(_x), "left");
  model.applyBC(BC::Dirichlet::FlagOnly(_x), "right");
  auto disp_increment = strain_rate * (posx_max - posx_min) / 2. * time_step;

  model.assembleInternalForces();

  for (UInt s = 1; s <= max_steps; ++s) {
    model.checkCohesiveStress();
    model.solveStep();

    auto nb_cohesive_elements = mesh.getNbElement(_cohesive_1d_2);

    if (s % 10 == 0) {
      std::cout << "passing step " << s << "/" << max_steps
                << ", number of cohesive elemets:" << nb_cohesive_elements
                << std::endl;
    }

    /// update external work and boundary conditions
    model.applyBC(BC::Dirichlet::IncrementValue(-disp_increment, _x), "left");
    model.applyBC(BC::Dirichlet::IncrementValue(disp_increment, _x), "right");
  }

  auto Ed = model.getEnergy("dissipated");
  auto Edt = 100. * 3.;

  std::cout << Ed << " " << Edt << std::endl;
  if (std::abs(Ed - Edt) > 0.001 || std::isnan(Ed)) {
    std::cout << "The dissipated energy is incorrect" << std::endl;
    finalize();
    return EXIT_FAILURE;
  }

  finalize();
  return EXIT_SUCCESS;
}
