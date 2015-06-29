/**
 * @file   test_solid_mechanics_model_igfem.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  test the solidmechancis model for IGFEM analysis
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model_igfem.hh"
#include "aka_common.hh"
// #include "mesh_segment_intersector.hh"
// #include "mesh_sphere_intersector.hh"
// #include "geom_helper_functions.hh"
// #include "mesh_geom_common.hh"
/* -------------------------------------------------------------------------- */
#include <cstdlib>
#include <fstream>
#include <iostream>
/* -------------------------------------------------------------------------- */
using namespace akantu;

int main(int argc, char *argv[]) {
  initialize("material.dat", argc, argv);
  UInt spatial_dimension = 2;
  /// mesh creation
  Mesh mesh(spatial_dimension);
  mesh.read("patch_test.msh");
  /// model creation
  SolidMechanicsModelIGFEM model(mesh);
  model.initFull();

  /// add fields that should be dumped
  model.setBaseName("IGFEM_test");
  model.addDumpField("material_index");
  model.dump();
  finalize();
  return EXIT_SUCCESS;
}
