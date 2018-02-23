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
#include "global_ids_updater.hh"
#include "mesh_accessor.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
void CohesiveElementInserter::initParallel(ElementSynchronizer & synchronizer) {
  global_ids_updater = std::make_unique<GlobalIdsUpdater>(mesh, synchronizer);
}

/* -------------------------------------------------------------------------- */
void CohesiveElementInserter::updateNodesType(Mesh & mesh,
                                              NewNodesEvent & node_event) {
  AKANTU_DEBUG_IN();

  auto & cohesive_node_event =
      dynamic_cast<CohesiveNewNodesEvent &>(node_event);
  auto & new_nodes = cohesive_node_event.getList();
  auto & old_nodes = cohesive_node_event.getOldNodesList();

  UInt local_nb_new_nodes = new_nodes.size();

  MeshAccessor mesh_accessor(mesh);
  auto & nodes_type = mesh_accessor.getNodesType();
  UInt nb_old_nodes = nodes_type.size();
  nodes_type.resize(nb_old_nodes + local_nb_new_nodes);

  for (auto && data : zip(old_nodes, new_nodes)) {
    UInt old_node, new_node;
    std::tie(old_node, new_node) = data;
    nodes_type(new_node) = nodes_type(old_node);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
UInt CohesiveElementInserter::updateGlobalIDs(NewNodesEvent & node_event) {
  AKANTU_DEBUG_IN();

  Array<UInt> & doubled_nodes = node_event.getList();

  UInt total_nb_new_nodes =
      global_ids_updater->updateGlobalIDsLocally(doubled_nodes.size());

  AKANTU_DEBUG_OUT();
  return total_nb_new_nodes;
}

void CohesiveElementInserter::synchronizeGlobalIDs(
    NewNodesEvent & /*node_event*/) {
  AKANTU_DEBUG_IN();

  global_ids_updater->synchronizeGlobalIDs();

  AKANTU_DEBUG_OUT();
}

} // namespace akantu
