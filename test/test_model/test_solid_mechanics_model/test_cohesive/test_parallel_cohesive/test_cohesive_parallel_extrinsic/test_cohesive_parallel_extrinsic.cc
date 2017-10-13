/**
 * @file   test_cohesive_parallel_extrinsic.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 *
 * @brief  parallel test for cohesive elements
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model_cohesive.hh"
#include "dumper_paraview.hh"
/* -------------------------------------------------------------------------- */
using namespace akantu;

int main(int argc, char *argv[]) {
  initialize("material.dat", argc, argv);

  const UInt max_steps = 500;

  UInt spatial_dimension = 2;
  Mesh mesh(spatial_dimension);

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  akantu::MeshPartition * partition = NULL;
  if(prank == 0) {
    // Read the mesh
    mesh.read("mesh.msh");

    /// partition the mesh
    partition = new MeshPartitionScotch(mesh, spatial_dimension);
    //    debug::setDebugLevel(dblDump);
    partition->partitionate(psize);
    //    debug::setDebugLevel(dblWarning);
  }

  SolidMechanicsModelCohesive model(mesh);

  model.initParallel(partition, NULL, true);

  // debug::setDebugLevel(dblDump);
  // std::cout << mesh << std::endl;
  // debug::setDebugLevel(dblWarning);


  model.initFull(SolidMechanicsModelCohesiveOptions(_explicit_lumped_mass, true));
  model.limitInsertion(_y, -0.30, -0.20);
  model.updateAutomaticInsertion();

  // debug::setDebugLevel(dblDump);
  // std::cout << mesh_facets << std::endl;
  // debug::setDebugLevel(dblWarning);

  Real time_step = model.getStableTimeStep()*0.1;
  model.setTimeStep(time_step);
  std::cout << "Time step: " << time_step << std::endl;

  model.assembleMassLumped();


  Array<Real> & position = mesh.getNodes();
  Array<Real> & velocity = model.getVelocity();
  Array<bool> & boundary = model.getBlockedDOFs();
  Array<Real> & displacement = model.getDisplacement();
  //  const Array<Real> & residual = model.getResidual();

  UInt nb_nodes = mesh.getNbNodes();

  /// boundary conditions
  for (UInt n = 0; n < nb_nodes; ++n) {
    if (position(n, 1) > 0.99 || position(n, 1) < -0.99)
      boundary(n, 1) = true;

    if (position(n, 0) > 0.99 || position(n, 0) < -0.99)
      boundary(n, 0) = true;
  }

  /// initial conditions
  Real loading_rate = 0.5;
  Real disp_update = loading_rate * time_step;
  for (UInt n = 0; n < nb_nodes; ++n) {
    velocity(n, 1) = loading_rate * position(n, 1);
  }

  model.synchronizeBoundaries();
  model.updateResidual();

  model.setBaseName("extrinsic_parallel");
  model.addDumpFieldVector("displacement");
  model.addDumpField("velocity"    );
  model.addDumpField("acceleration");
  model.addDumpField("residual"    );
  model.addDumpField("stress");
  model.addDumpField("grad_u");
  model.addDumpField("partitions");
  //  model.getDumper().getDumper().setMode(iohelper::BASE64);
  model.dump();

  model.setBaseNameToDumper("cohesive elements",
			    "extrinsic_parallel_cohesive_elements");
  model.addDumpFieldVectorToDumper("cohesive elements", "displacement");
  model.addDumpFieldToDumper("cohesive elements", "damage");
  model.dump("cohesive elements");

  /// Main loop
  for (UInt s = 1; s <= max_steps; ++s) {

    /// update displacement on extreme nodes
    for (UInt n = 0; n < nb_nodes; ++n) {
      if (position(n, 1) > 0.99 || position(n, 1) < -0.99)
	displacement(n, 1) += disp_update * position(n, 1);
    }

    model.checkCohesiveStress();
    model.solveStep();

    // model.dump();
    if(s % 10 == 0) {
      if(prank == 0) std::cout << "passing step " << s << "/" << max_steps << std::endl;
    }

    // // update displacement
    // for (UInt n = 0; n < nb_nodes; ++n) {
    //   if (position(n, 1) + displacement(n, 1) > 0) {
    // 	displacement(n, 0) -= 0.01;
    //   }
    // }

    //    Real Ed = dynamic_cast<MaterialCohesive&> (model.getMaterial(1)).getDissipatedEnergy();
    //    Real Er = dynamic_cast<MaterialCohesive&> (model.getMaterial(1)).getReversibleEnergy();

    // edis << s << " "
    // 	 << Ed << std::endl;

    // erev << s << " "
    // 	 << Er << std::endl;

  }

  model.dump();
  model.dump("cohesive elements");

  // edis.close();
  // erev.close();

  Real Ed = model.getEnergy("dissipated");

  Real Edt = 200 * sqrt(2);


  if(prank == 0) {
    std::cout << Ed << " " << Edt << std::endl;

    if (Ed < Edt * 0.999 || Ed > Edt * 1.001 || std::isnan(Ed)) {
      std::cout << "The dissipated energy is incorrect" << std::endl;
      finalize();
      return EXIT_FAILURE;
    }
  }

  finalize();
  return EXIT_SUCCESS;
}
