/**
 * @file   test_solid_mechanics_model_bar_traction2d_mass_not_lumped.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Tue Dec 20 15:37:09 2011
 *
 * @brief  test of the class SolidMechanicsModel
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include <limits>
#include <fstream>

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */
using namespace akantu;

akantu::ElementType type = akantu::_triangle_3;

akantu::SolidMechanicsModel * model;
akantu::UInt spatial_dimension = 2;
akantu::UInt nb_nodes;
akantu::UInt nb_element;

akantu::Array<akantu::Real> * lumped;

int main(int argc, char *argv[])
{
  akantu::initialize(argc, argv);
  akantu::UInt max_steps = 5000;
  akantu::Real time_factor = 0.8;

  //  akantu::Real epot, ekin;

  akantu::Mesh mesh(spatial_dimension);
  akantu::MeshIOMSH mesh_io;
  mesh_io.read("bar1.msh", mesh);

  model = new akantu::SolidMechanicsModel(mesh);

  nb_nodes = model->getFEEngine().getMesh().getNbNodes();
  nb_element = model->getFEEngine().getMesh().getNbElement(type);

  lumped = new akantu::Array<akantu::Real>(nb_nodes, spatial_dimension);

  /// model initialization
  model->initFull("material.dat");
  std::cout << model->getMaterial(0) << std::endl;

  model->initMaterials();
  model->initSolver();
  model->assembleMass();

  model->getMassMatrix().lump(*lumped);

  /// boundary conditions
  akantu::Real eps = 1e-16;
  const akantu::Array<akantu::Real> & pos = mesh.getNodes();
  akantu::Array<akantu::Real> & disp = model->getDisplacement();
  akantu::Array<bool> & boun = model->getBlockedDOFs();

  for (akantu::UInt i = 0; i < nb_nodes; ++i) {
    if(pos(i, 0) >= 9.) disp(i, 0) = (pos(i, 0) - 9) / 100.;
    if(pos(i) <= eps)   boun(i, 0) = true;
    if(pos(i, 1) <= eps || pos(i, 1) >= 1 - eps ) boun(i, 1) = true;
  }

  /// set the time step
  akantu::Real time_step = model->getStableTimeStep() * time_factor;
  std::cout << "Time Step = " << time_step << "s" << std::endl;
  model->setTimeStep(time_step);

  model->updateResidual();
  model->initialAcceleration();

  model->setBaseName("bar2d_mass_not_lumped");
  model->addDumpField("displacement");
  model->addDumpField("velocity"    );
  model->addDumpField("acceleration");
  model->addDumpField("force"       );
  model->addDumpField("residual"    );
  model->dump();

  std::ofstream energy;
  energy.open("energy_bar_2d_not_lumped.csv");
  energy << "id,rtime,epot,ekin,tot" << std::endl;

  for(akantu::UInt s = 1; s <= max_steps; ++s) {

    model->explicitPred();
    model->updateResidual();
    model->updateAcceleration();
    model->explicitCorr();

    akantu::Real epot = model->getPotentialEnergy();
    akantu::Real ekin = model->getKineticEnergy();

    energy << s << "," << (s-1)*time_step << "," << epot << "," << ekin << "," << epot + ekin
	   << std::endl;

    if(s % 1 == 0) model->dump();
    if(s % 100 == 0) std::cout << "passing step " << s << "/" << max_steps << std::endl;
  }

  energy.close();

  delete model;

  akantu::finalize();

  return EXIT_SUCCESS;
}
