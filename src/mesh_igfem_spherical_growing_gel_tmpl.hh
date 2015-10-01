/**
 * @file   mesh_igfem_spherical_growing_gel.cc
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Mon Sep 21 12:42:11 2015
 *
 * @brief  Implementation of the functions of MeshIgfemSphericalGrowingGel
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
#include "mesh_utils.hh"
#include <algorithm>

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template <UInt dim>
void MeshIgfemSphericalGrowingGel<dim>::buildSegmentConnectivityToNodeType() {
  Mesh mesh_facets(mesh.initMeshFacets());
  MeshUtils::buildSegmentToNodeType(mesh, mesh_facets, synchronizer);

  // only the ghost elements are considered
  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType ghost_type = (GhostType) g;

    Mesh::type_iterator it  = mesh_facets.firstType(1, ghost_type);
    Mesh::type_iterator end = mesh_facets.lastType(1, ghost_type);
    for(; it != end; ++it) {
      ElementType type = *it;

      const Array<Int> & segment_to_nodetype
	= mesh_facets.getData<Int>("segment_to_nodetype", type, ghost_type);

      const Array<UInt> & segment_connectivity
	= mesh_facets.getConnectivity(type, ghost_type);

      // looping over all the segments
      Array<UInt>::const_vector_iterator conn_it
	= segment_connectivity.begin(segment_connectivity.getNbComponent());
      Array<UInt>::const_vector_iterator conn_end
	= segment_connectivity.end(segment_connectivity.getNbComponent());
      UInt seg_index = 0;

      UInt connectivity_vals[2];
      for (; conn_it != conn_end; ++conn_it, ++seg_index) {
	connectivity_vals[0] = (*conn_it)(0);
	connectivity_vals[1] = (*conn_it)(1);

	if (!(mesh.isLocalNode(connectivity_vals[0]) &&
	      mesh.isLocalNode(connectivity_vals[1])) &&
	    !(mesh.isPureGhostNode(connectivity_vals[0]) &&
	      mesh.isPureGhostNode(connectivity_vals[1])) ) {
	  std::sort(connectivity_vals, connectivity_vals+2);

	  segment_conn_to_node_type[std::make_pair(connectivity_vals[0],
						   connectivity_vals[1])]
	    = segment_to_nodetype(seg_index);
	}
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void MeshIgfemSphericalGrowingGel<dim>::updateNodeType(const Array<UInt> & nodes_list,
						       const ElementTypeMapUInt & new_node_per_elem,
						       ElementType type) {
  Array<Int> & nodes_type = mesh.getNodesType();
  UInt old_nodes = nodes_type.getSize();

  // exit this function if the simulation in run in serial
  if (old_nodes == 0) return;

  UInt new_nodes = nodes_list.getSize();
  nodes_type.resize(old_nodes + new_nodes);
  Array<bool> checked_node(new_nodes, 1, false);

  UInt nb_prim_by_el = 0;
  if( (type == _triangle_3) ||
      (type == _igfem_triangle_4) ||
      (type == _igfem_triangle_5) ){
    nb_prim_by_el = 3;
  } else {
    AKANTU_DEBUG_ERROR("Not ready for mesh type " << type);
  }

  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType ghost_type = (GhostType) g;
    if (!this->mesh.getConnectivities().exists(type, ghost_type)) continue;

    // determine the node type for the local, master and slave nodes
    const Array<UInt> & connectivity = mesh.getConnectivity(type, ghost_type);
    const Array<UInt> & new_node_per_elem_array = new_node_per_elem(type, ghost_type);

    for (UInt el = 0; el < new_node_per_elem_array.getSize(); ++el) {
      UInt nb_nodes = new_node_per_elem_array(el, 0);

      for (UInt n = 0; n < nb_nodes; ++n) {
	UInt node_index = new_node_per_elem_array(el, 2*n+1);
	if (checked_node(node_index - old_nodes)) continue;

	// get the elemental connectivity of the segment associated to the node
	UInt segment_index = new_node_per_elem_array(el, 2*n+2) - 1;

	UInt extreme_nodes[2];
	extreme_nodes[0] = segment_index;
	extreme_nodes[1] = segment_index + 1;
	if (extreme_nodes[1] == nb_prim_by_el) extreme_nodes[1] = 0;

	// get the connectivity of the segment
	extreme_nodes[0] = connectivity(el, extreme_nodes[0]);
	extreme_nodes[1] = connectivity(el, extreme_nodes[1]);


	// if both extreme nodes are local, then also the new node is local
	if (mesh.isLocalNode(extreme_nodes[0]) &&
	    mesh.isLocalNode(extreme_nodes[1]))
	  nodes_type(node_index) = -1;
	// if both extreme nodes are pure ghost, then also the new node is pure ghost
	else if (mesh.isPureGhostNode(extreme_nodes[0]) &&
		 mesh.isPureGhostNode(extreme_nodes[1]))
	  nodes_type(node_index) = -3;
	// otherwise use the value stored in the map
	else {
	  std::sort(extreme_nodes, extreme_nodes+2);

	  AKANTU_DEBUG_ASSERT(segment_conn_to_node_type
			      .find(std::make_pair(extreme_nodes[0],
						   extreme_nodes[1]))
			      != segment_conn_to_node_type.end(),
			      "This segment should be present in the map at this point");

	  nodes_type(node_index)
	    = segment_conn_to_node_type[std::make_pair(extreme_nodes[0],
						       extreme_nodes[1])];
	}

	checked_node(node_index - old_nodes) = true;
      }
    }
  }

  AKANTU_DEBUG_ASSERT(std::accumulate(checked_node.begin(), checked_node.end(), 0)
		      == checked_node.getSize(),
		      "Not all new nodes were assigned a node type");
}

__END_AKANTU__
