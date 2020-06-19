/**
 * @file   new_local_material.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Aug 06 2015
 * @date last modification: Mon Jan 18 2016
 *
 * @brief  test of the class SolidMechanicsModel
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
#include "local_material_damage.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */
using namespace akantu;

#define bar_length 10.
#define bar_height 4.
akantu::Real eps = 1e-10;

int main(int argc, char * argv[]) {
  akantu::initialize("material.dat", argc, argv);

  UInt max_steps = 10000;
  Real epot, ekin;

  const UInt spatial_dimension = 2;

  Mesh mesh(spatial_dimension);
  mesh.read("barre_trou.msh");

  /// model creation
  SolidMechanicsModel model(mesh);

  /// model initialization
  model.initFull(_analysis_method = _explicit_lumped_mass);

  std::cout << model.getMaterial(0) << std::endl;

  Real time_step = model.getStableTimeStep();
  model.setTimeStep(time_step / 10.);

  /// Dirichlet boundary conditions
  model.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "Fixed_x");
  model.applyBC(BC::Dirichlet::FixedValue(0.0, _y), "Fixed_y");

  // Neumann boundary condition
  Matrix<Real> stress(2, 2);
  stress.eye(3e2);
  model.applyBC(BC::Neumann::FromStress(stress), "Traction");

  model.setBaseName("local_material");
  model.addDumpField("displacement");
  model.addDumpField("velocity");
  model.addDumpField("acceleration");
  model.addDumpField("external_force");
  model.addDumpField("internal_force");
  model.addDumpField("grad_u");
  model.addDumpField("stress");
  model.addDumpField("damage");
  model.dump();

  for (UInt s = 0; s < max_steps; ++s) {
    model.solveStep();

    epot = model.getEnergy("potential");
    ekin = model.getEnergy("kinetic");

    if (s % 100 == 0)
      std::cout << s << " " << epot << " " << ekin << " " << epot + ekin
                << std::endl;

    if (s % 1000 == 0)
      model.dump();
  }

  akantu::finalize();
  return EXIT_SUCCESS;
}
