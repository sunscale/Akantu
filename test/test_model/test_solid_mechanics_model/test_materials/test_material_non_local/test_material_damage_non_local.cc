/**
 * @file   test_material_damage_non_local.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Sep 13 2011
 * @date last modification: Thu Apr 03 2014
 *
 * @brief  test for non-local damage materials
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
#include <iostream>

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char *argv[])
{
  debug::setDebugLevel(dblWarning);

  akantu::initialize("material_damage_non_local.dat", argc, argv);
  UInt max_steps = 40000;

  const UInt spatial_dimension = 2;
  Mesh mesh(spatial_dimension);
  mesh.read("mesh.msh");
  mesh.createGroupsFromMeshData<std::string>("physical_names");

  SolidMechanicsModel model(mesh);

  /// model initialization

  model.initFull();

  Real time_step = model.getStableTimeStep();
  model.setTimeStep(time_step/10.);

  std::cout << model << std::endl;

  model.applyBC(BC::Dirichlet::FixedValue(0.0), "Fixed");

  // Boundary condition (Neumann)
  Matrix<Real> stress(2,2);
  stress.eye(3e6);
  model.applyBC(BC::Neumann::FromHigherDim(stress), "Traction");

  model.setBaseName("damage_non_local");
  model.addDumpField("displacement");
  model.addDumpField("mass"        );
  model.addDumpField("velocity"    );
  model.addDumpField("acceleration");
  model.addDumpField("force"       );
  model.addDumpField("residual"    );
  model.addDumpField("damage"      );
  model.addDumpField("stress"      );
  model.addDumpField("strain"      );
  model.dump();

  for(UInt s = 0; s < max_steps; ++s) {
    model.explicitPred();

    model.updateResidual();
    model.updateAcceleration();
    model.explicitCorr();

    if(s % 100 == 0) std::cout << "Step " << s+1 << "/" << max_steps <<std::endl;

    if(s % 100 == 0) model.dump();
  }

  akantu::finalize();
  return EXIT_SUCCESS;
}
