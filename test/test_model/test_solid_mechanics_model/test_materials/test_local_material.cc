/**
 * @file   test_local_material.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 *
 * @date creation: Fri Nov 26 2010
 * @date last modification: Fri Sep 19 2014
 *
 * @brief  test of the class SolidMechanicsModel
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "local_material_damage.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char *argv[])
{
  akantu::initialize("material.dat", argc, argv);
  UInt max_steps = 200;
  Real epot, ekin;

  const UInt spatial_dimension = 2;
  Mesh mesh(spatial_dimension);
  mesh.read("barre_trou.msh");
  mesh.createGroupsFromMeshData<std::string>("physical_names");

  SolidMechanicsModel model(mesh);

  /// model initialization
  model.initFull(SolidMechanicsModelOptions(_explicit_lumped_mass, true));
  model.registerNewCustomMaterials<LocalMaterialDamage>("local_damage");
  model.initMaterials();

  Real time_step = model.getStableTimeStep();
  model.setTimeStep(time_step/10.);

  model.assembleMassLumped();

  std::cout << model << std::endl;

  /// Dirichlet boundary conditions
  model.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "Fixed");
  model.applyBC(BC::Dirichlet::FixedValue(0.0, _y), "Fixed");

  // Boundary condition (Neumann)
  Matrix<Real> stress(2,2);
  stress.eye(3e6);
  model.applyBC(BC::Neumann::FromHigherDim(stress), "Traction");

  for(UInt s = 0; s < max_steps; ++s) {
    model.solveStep();

    epot = model.getPotentialEnergy();
    ekin = model.getKineticEnergy();

    std::cout << s << " " << epot << " " << ekin << " " << epot + ekin
	      << std::endl;
  }

  akantu::finalize();
  return EXIT_SUCCESS;
}
