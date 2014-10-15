/**
 * @file   solid_mechanics_model_cohesive_parallel.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date   Fri Jun 14 16:59:55 2013
 *
 * @brief  Functions for parallel cohesive elements
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
#include "solid_mechanics_model_cohesive.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::initParallel(MeshPartition * partition,
                                               DataAccessor * data_accessor,
                                               bool extrinsic) {
  AKANTU_DEBUG_IN();

  SolidMechanicsModel::initParallel(partition, data_accessor);

  /// create the distributed synchronizer for cohesive elements
  cohesive_distributed_synchronizer =
    new DistributedSynchronizer(mesh, "cohesive_distributed_synchronizer");

  synch_registry->registerSynchronizer(*cohesive_distributed_synchronizer,
  				       _gst_material_id);

  synch_registry->registerSynchronizer(*cohesive_distributed_synchronizer,
  				       _gst_smm_stress);

  synch_registry->registerSynchronizer(*cohesive_distributed_synchronizer,
  				       _gst_smm_boundary);

  synch_parallel->filterElementsByKind(cohesive_distributed_synchronizer,
				       _ek_cohesive);

  inserter = new CohesiveElementInserter(mesh, extrinsic, synch_parallel,
					 id+":cohesive_element_inserter");
  Mesh & mesh_facets = inserter->getMeshFacets();

  facet_synchronizer =
    FacetSynchronizer::createFacetSynchronizer(*synch_parallel,
					       mesh_facets);

  synch_registry->registerSynchronizer(*facet_synchronizer, _gst_smmc_facets);
  synch_registry->registerSynchronizer(*facet_synchronizer, _gst_smmc_facets_conn);

  synchronizeGhostFacetsConnectivity();

  /// create the facet synchronizer for extrinsic simulations
  if (extrinsic) {
    facet_stress_synchronizer =
      FacetStressSynchronizer::createFacetStressSynchronizer(*facet_synchronizer,
							     mesh_facets);

    synch_registry->registerSynchronizer(*facet_stress_synchronizer,
					 _gst_smmc_facets_stress);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::synchronizeGhostFacetsConnectivity() {
  AKANTU_DEBUG_IN();

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();

  if (psize > 1) {

    /// get global connectivity for not ghost facets
    global_connectivity = new ElementTypeMapArray<UInt>("global_connectivity", id);

    Mesh & mesh_facets = inserter->getMeshFacets();

    mesh_facets.initElementTypeMapArray(*global_connectivity, 1,
				       spatial_dimension - 1, true,
				       _ek_regular, true);

    mesh_facets.getGlobalConnectivity(*global_connectivity,
				      spatial_dimension - 1, _not_ghost);

    /// communicate
    synch_registry->synchronize(_gst_smmc_facets_conn);

    /// flip facets
    MeshUtils::flipFacets(mesh_facets, *global_connectivity, _ghost);

    delete global_connectivity;
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::updateCohesiveSynchronizers() {
  /// update synchronizers if needed

  if (facet_synchronizer != NULL)
    facet_synchronizer->updateDistributedSynchronizer(*cohesive_distributed_synchronizer,
						      *this,
                                                      mesh);

  if (facet_stress_synchronizer != NULL) {
    const ElementTypeMapArray<UInt> & prank_to_element
      = synch_parallel->getPrankToElement();

    facet_stress_synchronizer->updateFacetStressSynchronizer(*inserter,
							     prank_to_element,
							     *this);
  }
}

__END_AKANTU__
