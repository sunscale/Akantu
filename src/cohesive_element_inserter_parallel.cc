/**
 * @file   cohesive_element_inserter_parallel.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 *
 * @brief  Parallel functions for the cohesive element inserter
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "cohesive_element_inserter.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
void CohesiveElementInserter::initParallel(FacetSynchronizer * facet_synchronizer,
					   DistributedSynchronizer * distributed_synchronizer) {
  AKANTU_DEBUG_IN();

  this->facet_synchronizer = facet_synchronizer;

  global_ids_updater = new GlobalIdsUpdater(mesh, distributed_synchronizer);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void CohesiveElementInserter::updateNodesType(Mesh & mesh,
					      NewNodesEvent & node_event) {
  AKANTU_DEBUG_IN();

  Array<UInt> & doubled_nodes = node_event.getList();
  UInt local_nb_new_nodes = doubled_nodes.getSize();
  
  Array<Int> & nodes_type = mesh.getNodesType();
  UInt nb_old_nodes = nodes_type.getSize();
  nodes_type.resize(nb_old_nodes + local_nb_new_nodes);

  for (UInt n = 0; n < local_nb_new_nodes; ++n) {
    UInt old_node = doubled_nodes(n, 0);
    UInt new_node = doubled_nodes(n, 1);
    nodes_type(new_node) = nodes_type(old_node);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
UInt CohesiveElementInserter::updateGlobalIDs(NewNodesEvent & node_event) {
  AKANTU_DEBUG_IN();

  Array<UInt> & doubled_nodes = node_event.getList();

  UInt total_nb_new_nodes
    = global_ids_updater->updateGlobalIDsLocally(doubled_nodes.getSize());

  AKANTU_DEBUG_OUT();
  return total_nb_new_nodes;
}

void CohesiveElementInserter::synchronizeGlobalIDs(NewNodesEvent & node_event) {
  AKANTU_DEBUG_IN();

  global_ids_updater->synchronizeGlobalIDs();

  AKANTU_DEBUG_OUT();
}



__END_AKANTU__
