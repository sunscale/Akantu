/**
 * @file   test_facet_synchronizer.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 *
 * @brief  Facet synchronizer test
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "facet_synchronizer.hh"
#include "mesh_utils.hh"
#include "synchronizer_registry.hh"
#include "test_data_accessor.hh"
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char * argv[]) {
  akantu::initialize(argc, argv);

  UInt spatial_dimension = 3;
  Mesh mesh(spatial_dimension);

  const auto & comm = Communicator::getStaticCommunicator();
  Int prank = comm.whoAmI();

  /// partition the mesh
  if (prank == 0) {
    mesh.read("facet.msh");
  }

  mesh.distribute();

  auto & synchronizer = mesh.getElementSynchronizer();

  /// create facets
  Mesh & mesh_facets = mesh.initMeshFacets("mesh_facets");

  /// setup facet communications
  auto & facet_synchronizer = mesh_facets.getElementSynchronizer();

  // AKANTU_DEBUG_INFO("Creating TestAccessor");
  // TestAccessor test_accessor(mesh);
  // SynchronizerRegistry synch_registry(test_accessor);

  // synch_registry.registerSynchronizer(facet_synchronizer, _gst_test);

  // /// synchronize facets and check results
  // AKANTU_DEBUG_INFO("Synchronizing tag");
  // synch_registry.synchronize(_gst_test);

  return 0;
}
