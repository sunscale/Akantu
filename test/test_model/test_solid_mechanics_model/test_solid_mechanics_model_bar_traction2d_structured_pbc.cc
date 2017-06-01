/**
 * @file   test_solid_mechanics_model_bar_traction2d_structured_pbc.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Dec 22 2011
 * @date last modification: Thu Aug 06 2015
 *
 * @brief  test of pbc class for SolidMechanicsModel
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
#include <limits>
#include <fstream>
#include <iostream>

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[])
{
  akantu::debug::setDebugLevel(akantu::dblWarning);
  akantu::initialize("material.dat", argc, argv);
  akantu::UInt spatial_dimension = 2;
  akantu::UInt max_steps = 1000;
  akantu::Real time_factor = 0.2;

  akantu::Real epot, ekin;

  akantu::Mesh mesh(spatial_dimension);
  mesh.read("bar_structured1.msh");

  akantu::SolidMechanicsModel model(mesh);

  /// model initialization
  model.initFull();

  std::cout << model.getMaterial(0) << std::endl;

  model.setPBC(1,0,0);
  model.initPBC();
  model.assembleMassLumped();

  model.setBaseName("bar2d_structured_pbc");
  model.addDumpField("displacement");
  model.addDumpField("mass"        );
  model.addDumpField("velocity"    );
  model.addDumpField("acceleration");
  model.addDumpFieldVector("external_force");
  model.addDumpFieldVector("internal_force");
  model.addDumpField("stress"      );
  model.addDumpField("strain"      );
  model.dump();

  akantu::UInt nb_nodes = mesh.getNbNodes();

  /// boundary conditions
  mesh.computeBoundingBox();
  akantu::Real eps = 1e-16;
  akantu::Real signal_start = 0.6*mesh.getUpperBounds()(0);
  akantu::Real signal_end = 0.7*mesh.getUpperBounds()(0);
  akantu::Real delta_d = signal_end - signal_start;
  akantu::Real signal = 1.;
  const akantu::Array<akantu::Real> & coords = model.getFEEngine().getMesh().getNodes();
  akantu::Array<akantu::Real> & disp = model.getDisplacement();
  for (akantu::UInt i = 0; i < nb_nodes; ++i) {
    if(coords(i,0) >= signal_start && coords(i,0) <= signal_end) {
      if (coords(i,0) <= 0.5 * (signal_start + signal_end))
	disp(i,0) = (coords(i,0) - signal_start) * 2 * signal / delta_d;
      else
	disp(i,0) = (signal_end - coords(i,0)) * 2 * signal / delta_d;
    }

    if(coords(i,1) <= eps || coords(i,1) >= 1 - eps ) {
      model.getBlockedDOFs().storage()[spatial_dimension*i + 1] = true;
    }
  }

  akantu::Real time_step = model.getStableTimeStep() * time_factor;
  std::cout << "Time Step = " << time_step << "s" << std::endl;
  model.setTimeStep(time_step);

  std::ofstream energy;
  energy.open("energy_2d_pbc.csv");
  energy << "id,epot,ekin,tot" << std::endl;

  for(akantu::UInt s = 1; s <= max_steps; ++s) {
    model.explicitPred();
    model.updateResidual();
    model.updateAcceleration();
    model.explicitCorr();

    epot = model.getEnergy("potential");
    ekin = model.getEnergy("kinetic");

    if(s % 20 == 0) model.dump();

    std::cerr << "passing step " << s << "/" << max_steps << std::endl;
    energy << s << "," << epot << "," << ekin << "," << epot + ekin
	   << std::endl;
  }

  energy.close();

  return EXIT_SUCCESS;
}
