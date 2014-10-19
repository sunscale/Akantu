/**
 * @file   test_solid_mechanics_model_segment_parallel.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sun Sep 26 2010
 * @date last modification: Thu Jun 05 2014
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
#include <limits>
#include <fstream>

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"

/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
  akantu::UInt spatial_dimension = 1;
  akantu::UInt max_steps = 10000;
  akantu::Real time_factor = 0.2;

  akantu::initialize("material.dat", argc, argv);

  akantu::Mesh mesh(spatial_dimension);
  akantu::StaticCommunicator & comm = akantu::StaticCommunicator::getStaticCommunicator();
  akantu::Int psize = comm.getNbProc();
  akantu::Int prank = comm.whoAmI();

  akantu::debug::setDebugLevel(akantu::dblInfo);

  akantu::MeshPartition * partition = NULL;
  if(prank == 0) {
    mesh.read("segment.msh");
    partition = new akantu::MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);
  }


  akantu::SolidMechanicsModel model(mesh);
  model.initParallel(partition);

  /// model initialization
  model.initFull();

  std::cout << model.getMaterial(0) << std::endl;

  model.assembleMassLumped();

  model.setBaseName("segment_parallel");
  model.addDumpField("displacement");
  model.addDumpField("mass"        );
  model.addDumpField("velocity"    );
  model.addDumpField("acceleration");
  model.addDumpField("force"       );
  model.addDumpField("residual"    );
  model.addDumpField("stress"      );
  model.addDumpField("strain"      );

  /// boundary conditions
  for (akantu::UInt i = 0; i < mesh.getNbNodes(); ++i) {
    model.getDisplacement().storage()[spatial_dimension*i] = model.getFEEngine().getMesh().getNodes().storage()[i] / 100. ;

    if(model.getFEEngine().getMesh().getNodes().storage()[spatial_dimension*i] <= 1e-15)
      model.getBlockedDOFs().storage()[i] = true;
  }

  akantu::Real time_step = model.getStableTimeStep() * time_factor;
  std::cout << "Time Step = " << time_step << "s" << std::endl;
  model.setTimeStep(time_step);

  for(akantu::UInt s = 1; s <= max_steps; ++s) {
    model.explicitPred();

    model.updateResidual();
    model.updateAcceleration();
    model.explicitCorr();

    akantu::Real epot = model.getPotentialEnergy();
    akantu::Real ekin = model.getKineticEnergy();

    if(prank == 0) {
      std::cout << s << " " << epot << " " << ekin << " " << epot + ekin
		<< std::endl;
    }

    model.dump();
    if(s % 10 == 0) std::cerr << "passing step " << s << "/" << max_steps << std::endl;
  }

  akantu::finalize();

  return EXIT_SUCCESS;
}
