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
void CohesiveElementInserter::initParallel(FacetSynchronizer * facet_synchronizer) {
  AKANTU_DEBUG_IN();

  this->facet_synchronizer = facet_synchronizer;

  distributed_synchronizer
    = new DistributedSynchronizer(mesh, "inserter_synchronizer");

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
  UInt local_nb_new_nodes = doubled_nodes.getSize();

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int rank = comm.whoAmI();
  Int nb_proc = comm.getNbProc();

  /// resize global ids array
  Array<UInt> & nodes_global_ids = *mesh.nodes_global_ids;
  UInt nb_old_nodes = nodes_global_ids.getSize();

  nodes_global_ids.resize(nb_old_nodes + local_nb_new_nodes);

  /// compute amount of local or master doubled nodes
  Vector<UInt> local_master_nodes(nb_proc);

  for (UInt n = 0; n < local_nb_new_nodes; ++n) {
    UInt old_node = doubled_nodes(n, 0);
    if (mesh.isLocalOrMasterNode(old_node)) ++local_master_nodes(rank);
  }

  comm.allGather(local_master_nodes.storage(), 1);

  /// update global number of nodes
  UInt total_nb_new_nodes = std::accumulate(local_master_nodes.storage(),
					    local_master_nodes.storage() + nb_proc,
					    0);

  if (total_nb_new_nodes == 0) return 0;

  /// set global ids of local and master nodes
  UInt starting_index = std::accumulate(local_master_nodes.storage(),
  					local_master_nodes.storage() + rank,
  					mesh.nb_global_nodes);

  for (UInt n = 0; n < local_nb_new_nodes; ++n) {
    UInt new_node = doubled_nodes(n, 1);
    if (mesh.isLocalOrMasterNode(new_node)) {
      nodes_global_ids(new_node) = starting_index;
      ++starting_index;
    }
  }

  /// update distributed synchronizer
  distributed_synchronizer->reset();
  facet_synchronizer->updateDistributedSynchronizer(*distributed_synchronizer,
						    *this, mesh);
  distributed_synchronizer->computeAllBufferSizes(*this);

  /// communicate global ids
  distributed_synchronizer->asynchronousSynchronize(*this, _gst_ce_inserter);
  distributed_synchronizer->waitEndSynchronize(*this, _gst_ce_inserter);

  AKANTU_DEBUG_OUT();
  return total_nb_new_nodes;
}


__END_AKANTU__
