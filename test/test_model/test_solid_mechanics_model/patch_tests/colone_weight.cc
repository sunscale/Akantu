/**
 * @file   colone_weight.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Alodie Schneuwly <alodie.schneuwly@epfl.ch>
 *
 * @date creation: Mon Aug 09 2010
 * @date last modification: Thu Aug 06 2015
 *
 * @brief  column test
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


/* -------------------------------------------------------------------------- */
#include <limits>
#include <fstream>
#include <iostream>
/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char *argv[]) {
  // chose if you use hexahedron elements
  bool use_hexa = false;

  std::stringstream mesh_file;
  std::stringstream output;
  std::stringstream energy_file;
  akantu::ElementType type;
   UInt vel_damping_interval;

  if (use_hexa) {
    type = akantu::_hexahedron_8;
    mesh_file << "colone_hexa.msh";
    output << "paraview/test_weight_hexa";
    energy_file << "energy_hexa.csv";
    vel_damping_interval =4;
  }
  else {
    type = akantu::_tetrahedron_4;
    mesh_file << "colone_tetra.msh";
    output << "paraview/test_weight_tetra";
    energy_file << "energy_tetra.csv";
    vel_damping_interval = 8;
  }

  akantu::UInt spatial_dimension = 3;

  akantu::UInt max_steps = 2000;
  akantu::Real time_factor = 0.8;

  akantu::initialize("material_colone.dat", argc, argv);

  //  akantu::Real epot, ekin;
  akantu::Mesh mesh(spatial_dimension);
  mesh.read(mesh_file.str().c_str());

  akantu::SolidMechanicsModel model(mesh);

  akantu::UInt nb_nodes = mesh.getNbNodes();
  akantu::UInt nb_element = mesh.getNbElement(type);

  std::cout << "Nb nodes : " << nb_nodes << " - nb elements : " << nb_element << std::endl;

  /// model initialization
  model.initFull();

  std::cout << model.getMaterial(0) << std::endl;

  model.assembleMassLumped();

  /// boundary conditions
  const akantu::Array<Real> & position = model.getFEEngine().getMesh().getNodes();
  akantu::Array<bool> & boundary = model.getBlockedDOFs();
  akantu::Array<Real> & force = model.getForce();
  const akantu::Array<Real> & mass = model.getMass();

  akantu::Real z_min = position(0, 2);
  for (unsigned int i = 0; i < nb_nodes; ++i) {
    if(position(i, 2) < z_min)
      z_min = position(i, 2);
  }

  akantu::Real eps = 1e-13;
  for (akantu::UInt i = 0; i < nb_nodes; ++i) {
    if(fabs(position(i, 2) - z_min) <= eps)
      boundary(i,2) = true;
    else
      force(i,2) = -mass(i,0) * 9.81;
  }

  akantu::Real time_step = model.getStableTimeStep() * time_factor;
  std::cout << "Time Step = " << time_step << "s" << std::endl;
  model.setTimeStep(time_step);

  model.updateResidual();

  model.setBaseName("colonne_weight");
  model.addDumpField("displacement");
  model.addDumpField("mass"        );
  model.addDumpField("velocity"    );
  model.addDumpField("acceleration");
  model.addDumpField("external_force");
  model.addDumpField("internal_force");
  model.addDumpField("damage"      );
  model.addDumpField("stress"      );
  model.addDumpField("strain"      );
  model.dump();

  akantu::Array<Real> & velocity = model.getVelocity();

  std::ofstream energy;
  energy.open(energy_file.str().c_str());
  energy << "id,epot,ekin,tot" << std::endl;

  for(akantu::UInt s = 1; s <= max_steps; ++s) {

    model.explicitPred();
    model.updateResidual();
    model.updateAcceleration();
    model.explicitCorr();

    akantu::Real epot = model.getEnergy("potential");
    akantu::Real ekin = model.getEnergy("kinetic");
    energy << s << "," << epot << "," << ekin << "," << epot + ekin
	   << std::endl;

    if (s % vel_damping_interval == 0) {
      for (akantu::UInt i = 0; i < nb_nodes; ++i) {
	velocity(i, 0) *= 0.9;
	velocity(i, 1) *= 0.9;
	velocity(i, 2) *= 0.9;
      }
    }

    if(s % 1 == 0) model.dump();
    if(s % 10 == 0) std::cout << "passing step " << s << "/" << max_steps << std::endl;
  }

  akantu::finalize();

  return EXIT_SUCCESS;
}
