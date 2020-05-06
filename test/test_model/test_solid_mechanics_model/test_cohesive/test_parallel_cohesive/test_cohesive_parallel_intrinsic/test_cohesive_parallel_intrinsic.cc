/**
 * @file   test_cohesive_parallel_intrinsic.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Nov 05 2014
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  parallel test for intrinsic cohesive elements
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

  const UInt max_steps = 350;

  UInt spatial_dimension = 2;
  Mesh mesh(spatial_dimension);

  const auto & comm = Communicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  akantu::MeshPartition * partition = NULL;
  if (prank == 0) {
    // Read the mesh
    mesh.read("mesh.msh");

    // /// insert cohesive elements
    // CohesiveElementInserter inserter(mesh);
    // inserter.setLimit('x', -0.26, -0.24);
    // inserter.insertIntrinsicElements();

    /// partition the mesh
    partition = new MeshPartitionScotch(mesh, spatial_dimension);
    //    debug::setDebugLevel(dblDump);
    partition->partitionate(psize);
    //   debug::setDebugLevel(dblWarning);
  }

  SolidMechanicsModelCohesive model(mesh);

  model.initParallel(partition);

  model.initFull();

  model.limitInsertion(_x, -0.26, -0.24);
  model.insertIntrinsicElements();

  debug::setDebugLevel(dblDump);
  std::cout << mesh << std::endl;
  debug::setDebugLevel(dblWarning);

  Real time_step = model.getStableTimeStep() * 0.8;
  model.setTimeStep(time_step);
  //  std::cout << "Time step: " << time_step << std::endl;

  model.assembleMassLumped();

  Array<Real> & position = mesh.getNodes();
  Array<Real> & velocity = model.getVelocity();
  Array<bool> & boundary = model.getBlockedDOFs();
  //  Array<Real> & displacement = model.getDisplacement();
  //  const Array<Real> & residual = model.getResidual();

  UInt nb_nodes = mesh.getNbNodes();
  Real epsilon = std::numeric_limits<Real>::epsilon();

  for (UInt n = 0; n < nb_nodes; ++n) {
    if (std::abs(position(n, 0) - 1.) < epsilon)
      boundary(n, 0) = true;
  }

  model.synchronizeBoundaries();
  model.updateResidual();

  model.setBaseName("intrinsic_parallel");
  model.addDumpFieldVector("displacement");
  model.addDumpField("velocity");
  model.addDumpField("acceleration");
  model.addDumpField("residual");
  model.addDumpField("stress");
  model.addDumpField("strain");
  model.addDumpField("partitions");
  model.addDumpField("force");
  model.dump();

  model.setBaseNameToDumper("cohesive elements",
                            "cohesive_elements_parallel_intrinsic");
  model.addDumpFieldVectorToDumper("cohesive elements", "displacement");
  model.dump("cohesive elements");

  /// initial conditions
  Real loading_rate = .2;
  for (UInt n = 0; n < nb_nodes; ++n) {
    velocity(n, 0) = loading_rate * position(n, 0);
  }

  /// Main loop
  for (UInt s = 1; s <= max_steps; ++s) {

    model.solveStep();

    if (s % 20 == 0) {
      model.dump();
      model.dump("cohesive elements");
      if (prank == 0)
        std::cout << "passing step " << s << "/" << max_steps << std::endl;
    }

    // // update displacement
    // for (UInt n = 0; n < nb_nodes; ++n) {
    //   if (position(n, 1) + displacement(n, 1) > 0) {
    // 	displacement(n, 0) -= 0.01;
    //   }
    // }

    //    Real Ed = dynamic_cast<MaterialCohesive&>
    //    (model.getMaterial(1)).getDissipatedEnergy();
    //    Real Er = dynamic_cast<MaterialCohesive&>
    //    (model.getMaterial(1)).getReversibleEnergy();

    // edis << s << " "
    // 	 << Ed << std::endl;

    // erev << s << " "
    // 	 << Er << std::endl;
  }

  // edis.close();
  // erev.close();

  Real Ed = model.getEnergy("dissipated");

  Real Edt = 2 * sqrt(2);

  if (prank == 0) {
    std::cout << Ed << " " << Edt << std::endl;

    if (std::abs((Ed - Edt) / Edt) > 0.01 || std::isnan(Ed)) {
      std::cout << "The dissipated energy is incorrect" << std::endl;
      return EXIT_FAILURE;
    }
  }

  finalize();
  if (prank == 0)
    std::cout << "OK: Test passed!" << std::endl;
  return EXIT_SUCCESS;
}
