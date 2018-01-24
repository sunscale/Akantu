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

#include "test_data_accessor.hh"
#include "facet_synchronizer.hh"
#include "mesh_utils.hh"
#include "synchronizer_registry.hh"
#include <iostream>

/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char *argv[]) {
  akantu::initialize(argc, argv);

  UInt spatial_dimension = 3;
  Mesh mesh(spatial_dimension);

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();
  DistributedSynchronizer * dist = NULL;

  /// partition the mesh
  if(prank == 0) {
    mesh.read("facet.msh");
    MeshPartition * partition = new MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);
    dist = DistributedSynchronizer::createDistributedSynchronizerMesh(mesh, partition);
    delete partition;
  } else {
    dist = DistributedSynchronizer::createDistributedSynchronizerMesh(mesh, NULL);
  }

  /// create facets
  Mesh mesh_facets(mesh.initMeshFacets("mesh_facets"));
  MeshUtils::buildAllFacets(mesh, mesh_facets, 0, dist);

  /// compute barycenter for each facet
  ElementTypeMapArray<Real> barycenters("barycenters", "", 0);
  mesh_facets.initElementTypeMapArray(barycenters,
				     spatial_dimension,
				     spatial_dimension - 1);

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {
    GhostType ghost_type = *gt;
    Mesh::type_iterator it = mesh_facets.firstType(spatial_dimension - 1,
						   ghost_type);
    Mesh::type_iterator last_type = mesh_facets.lastType(spatial_dimension - 1,
							 ghost_type);

    for(; it != last_type; ++it) {
      UInt nb_element = mesh_facets.getNbElement(*it, ghost_type);
      Array<Real> & barycenter = barycenters(*it, ghost_type);
      barycenter.resize(nb_element);

      Array<Real>::iterator< Vector<Real> > bary_it
	= barycenter.begin(spatial_dimension);

      for (UInt elem = 0; elem < nb_element; ++elem, ++bary_it)
	mesh_facets.getBarycenter(elem, *it, bary_it->storage(), ghost_type);
    }
  }

  /// setup facet communications
  FacetSynchronizer * facet_synchronizer =
    FacetSynchronizer::createFacetSynchronizer(*dist, mesh_facets);

  AKANTU_DEBUG_INFO("Creating TestAccessor");
  TestAccessor test_accessor(mesh, barycenters);
  SynchronizerRegistry synch_registry(test_accessor);

  synch_registry.registerSynchronizer(*facet_synchronizer, _gst_test);

  /// synchronize facets and check results
  AKANTU_DEBUG_INFO("Synchronizing tag");
  synch_registry.synchronize(_gst_test);

  delete facet_synchronizer;
  delete dist;
  akantu::finalize();

  return EXIT_SUCCESS;
}
