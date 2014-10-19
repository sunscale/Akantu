/**
 * @file   distributed_synchronizer.cc
 *
 * @author Dana Christen <dana.christen@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Jun 16 2011
 * @date last modification: Fri Sep 05 2014
 *
 * @brief  implementation of a  communicator using a static_communicator for real
 * send/receive
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_common.hh"
#include "distributed_synchronizer.hh"
#include "static_communicator.hh"
#include "mesh_utils.hh"
#include "mesh_data.hh"
#include "element_group.hh"
/* -------------------------------------------------------------------------- */
#include <map>
#include <iostream>
#include <algorithm>

#if defined(AKANTU_DEBUG_TOOLS)
#  include "aka_debug_tools.hh"
#endif

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
DistributedSynchronizer::DistributedSynchronizer(Mesh & mesh,
						 SynchronizerID id,
						 MemoryID memory_id) :
  Synchronizer(id, memory_id),
  mesh(mesh),
  static_communicator(&StaticCommunicator::getStaticCommunicator()),
  prank_to_element("prank_to_element", id)
{
  AKANTU_DEBUG_IN();

  nb_proc = static_communicator->getNbProc();
  rank    = static_communicator->whoAmI();

  send_element = new Array<Element>[nb_proc];
  recv_element = new Array<Element>[nb_proc];

  mesh.registerEventHandler(*this);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
DistributedSynchronizer::~DistributedSynchronizer() {
  AKANTU_DEBUG_IN();

  for (UInt p = 0; p < nb_proc; ++p) {
    send_element[p].clear();
    recv_element[p].clear();
  }

  delete [] send_element;
  delete [] recv_element;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
DistributedSynchronizer * DistributedSynchronizer::
createDistributedSynchronizerMesh(Mesh & mesh,
				  const MeshPartition * partition,
				  UInt root,
				  SynchronizerID id,
				  MemoryID memory_id) {
  AKANTU_DEBUG_IN();

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  UInt nb_proc = comm.getNbProc();
  UInt my_rank = comm.whoAmI();

  DistributedSynchronizer & communicator = *(new DistributedSynchronizer(mesh, id, memory_id));

  if(nb_proc == 1) return &communicator;

  UInt * local_connectivity = NULL;
  UInt * local_partitions = NULL;
  Array<UInt> * old_nodes = mesh.getNodesGlobalIdsPointer();
  old_nodes->resize(0);
  Array<Real> * nodes = mesh.getNodesPointer();


  UInt spatial_dimension = nodes->getNbComponent();

  mesh.synchronizeGroupNames();

  /* ------------------------------------------------------------------------ */
  /*  Local (rank == root)                                                    */
  /* ------------------------------------------------------------------------ */
  if(my_rank == root) {
    AKANTU_DEBUG_ASSERT(partition->getNbPartition() == nb_proc,
			"The number of partition does not match the number of processors: " <<
			partition->getNbPartition() << " != " << nb_proc);

    /**
     * connectivity and communications scheme construction
     */
    Mesh::type_iterator it  = mesh.firstType(_all_dimensions,
					     _not_ghost,
					     _ek_not_defined);
    Mesh::type_iterator end = mesh.lastType(_all_dimensions,
					    _not_ghost,
					    _ek_not_defined);
    UInt count = 0;
    /* --- MAIN LOOP ON TYPES --- */
    for(; it != end; ++it) {
      ElementType type = *it;

      UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
      UInt nb_element = mesh.getNbElement(*it);
      UInt nb_local_element[nb_proc];
      UInt nb_ghost_element[nb_proc];
      UInt nb_element_to_send[nb_proc];

      memset(nb_local_element, 0, nb_proc*sizeof(UInt));
      memset(nb_ghost_element, 0, nb_proc*sizeof(UInt));
      memset(nb_element_to_send, 0, nb_proc*sizeof(UInt));

      /// \todo change this ugly way to avoid a problem if an element
      /// type is present in the mesh but not in the partitions
      const Array<UInt> * tmp_partition_num = NULL;
      try {
	tmp_partition_num = &partition->getPartition(type, _not_ghost);
      } catch(...) {
	continue;
      }
      const Array<UInt> & partition_num = *tmp_partition_num;

      const CSR<UInt> & ghost_partition = partition->getGhostPartitionCSR()(type, _not_ghost);

      /* -------------------------------------------------------------------- */
      /// constructing the reordering structures
      for (UInt el = 0; el < nb_element; ++el) {
	nb_local_element[partition_num(el)]++;
	for (CSR<UInt>::const_iterator part = ghost_partition.begin(el);
	     part != ghost_partition.end(el);
	     ++part) {
	  nb_ghost_element[*part]++;
	}
	nb_element_to_send[partition_num(el)] += ghost_partition.getNbCols(el) + 1;
      }

      /// allocating buffers
      UInt * buffers[nb_proc];
      UInt * buffers_tmp[nb_proc];
      for (UInt p = 0; p < nb_proc; ++p) {
	UInt size = nb_nodes_per_element * (nb_local_element[p] +
					    nb_ghost_element[p]);
	buffers[p] = new UInt[size];
	buffers_tmp[p] = buffers[p];
      }

      /// copying the local connectivity
      UInt * conn_val = mesh.getConnectivity(type, _not_ghost).storage();
      for (UInt el = 0; el < nb_element; ++el) {
	memcpy(buffers_tmp[partition_num(el)],
	       conn_val + el * nb_nodes_per_element,
	       nb_nodes_per_element * sizeof(UInt));
	buffers_tmp[partition_num(el)] += nb_nodes_per_element;
      }

      /// copying the connectivity of ghost element
      for (UInt el = 0; el < nb_element; ++el) {
	for (CSR<UInt>::const_iterator part = ghost_partition.begin(el);
	     part != ghost_partition.end(el);
	     ++part) {
	  UInt proc = *part;
	  memcpy(buffers_tmp[proc],
		 conn_val + el * nb_nodes_per_element,
		 nb_nodes_per_element * sizeof(UInt));
	  buffers_tmp[proc] += nb_nodes_per_element;
	}
      }


      /// tag info
      std::vector<std::string> tag_names;
      mesh.getMeshData().getTagNames(tag_names, type);

      UInt nb_tags = tag_names.size();

      /* -------->>>>-SIZE + CONNECTIVITY------------------------------------ */
      /// send all connectivity and ghost information to all processors
      std::vector<CommunicationRequest *> requests;
      for (UInt p = 0; p < nb_proc; ++p) {
	if(p != root) {
	  UInt size[5];
	  size[0] = (UInt) type;
	  size[1] = nb_local_element[p];
	  size[2] = nb_ghost_element[p];
	  size[3] = nb_element_to_send[p];
	  size[4] = nb_tags;
	  AKANTU_DEBUG_INFO("Sending connectivities informations to proc " << p << " TAG("<< Tag::genTag(my_rank, count, TAG_SIZES) <<")");
	  comm.send(size, 5, p, Tag::genTag(my_rank, count, TAG_SIZES));

	  AKANTU_DEBUG_INFO("Sending connectivities to proc " << p << " TAG("<< Tag::genTag(my_rank, count, TAG_CONNECTIVITY) <<")");
	  requests.push_back(comm.asyncSend(buffers[p],
					    nb_nodes_per_element * (nb_local_element[p] +
								    nb_ghost_element[p]),
					    p, Tag::genTag(my_rank, count, TAG_CONNECTIVITY)));
	} else {
	  local_connectivity = buffers[p];
	}
      }

      /// create the renumbered connectivity
      AKANTU_DEBUG_INFO("Renumbering local connectivities");
      MeshUtils::renumberMeshNodes(mesh,
				   local_connectivity,
				   nb_local_element[root],
				   nb_ghost_element[root],
				   type,
				   *old_nodes);

      comm.waitAll(requests);
      comm.freeCommunicationRequest(requests);
      requests.clear();

      for (UInt p = 0; p < nb_proc; ++p) {
	delete [] buffers[p];
      }

      /* -------------------------------------------------------------------- */
      for (UInt p = 0; p < nb_proc; ++p) {
	buffers[p] = new UInt[nb_ghost_element[p] + nb_element_to_send[p]];
	buffers_tmp[p] = buffers[p];
      }

      /// splitting the partition information to send them to processors
      UInt count_by_proc[nb_proc];
      memset(count_by_proc, 0, nb_proc*sizeof(UInt));
      for (UInt el = 0; el < nb_element; ++el) {
	*(buffers_tmp[partition_num(el)]++) = ghost_partition.getNbCols(el);
	UInt i(0);
	for (CSR<UInt>::const_iterator part = ghost_partition.begin(el);
	     part != ghost_partition.end(el);
	     ++part, ++i) {
	  *(buffers_tmp[partition_num(el)]++) = *part;
	}
      }

      for (UInt el = 0; el < nb_element; ++el) {
	UInt i(0);
	for (CSR<UInt>::const_iterator part = ghost_partition.begin(el);
	     part != ghost_partition.end(el);
	     ++part, ++i) {
	  *(buffers_tmp[*part]++) = partition_num(el);
	}
      }

      /* -------->>>>-PARTITIONS--------------------------------------------- */
      /// last data to compute the communication scheme
      for (UInt p = 0; p < nb_proc; ++p) {
	if(p != root) {
	  AKANTU_DEBUG_INFO("Sending partition informations to proc " << p << " TAG("<< Tag::genTag(my_rank, count, TAG_PARTITIONS) <<")");
	  requests.push_back(comm.asyncSend(buffers[p],
					     nb_element_to_send[p] + nb_ghost_element[p],
					    p, Tag::genTag(my_rank, count, TAG_PARTITIONS)));
	} else {
	  local_partitions = buffers[p];
	}
      }

      if(Mesh::getSpatialDimension(type) == mesh.getSpatialDimension()) {
	AKANTU_DEBUG_INFO("Creating communications scheme");
	communicator.fillCommunicationScheme(local_partitions,
					     nb_local_element[root],
					     nb_ghost_element[root],
					     type);
      }

      comm.waitAll(requests);
      comm.freeCommunicationRequest(requests);
      requests.clear();

      for (UInt p = 0; p < nb_proc; ++p) {
	delete [] buffers[p];
      }

      /* -------------------------------------------------------------------- */
      /// send  data assossiated to the mesh
      /* -------->>>>-TAGS--------------------------------------------------- */
      synchronizeTagsSend(communicator, root, mesh, nb_tags, type,
			  partition_num,
			  ghost_partition,
			  nb_local_element[root],
			  nb_ghost_element[root]);


      /* -------------------------------------------------------------------- */
      /// send  data assossiated to groups
      /* -------->>>>-GROUPS------------------------------------------------- */
      synchronizeElementGroups(communicator, root, mesh, type,
      			       partition_num,
      			       ghost_partition,
			       nb_element);

      ++count;
    }

    /* -------->>>>-SIZE----------------------------------------------------- */
    for (UInt p = 0; p < nb_proc; ++p) {
      if(p != root) {
	UInt size[5];
	size[0] = (UInt) _not_defined;
	size[1] = 0;
	size[2] = 0;
	size[3] = 0;
	size[4] = 0;
	AKANTU_DEBUG_INFO("Sending empty connectivities informations to proc " << p << " TAG("<< Tag::genTag(my_rank, count, TAG_SIZES) <<")");
	comm.send(size, 5, p, Tag::genTag(my_rank, count, TAG_SIZES));
      }
    }

    /* ---------------------------------------------------------------------- */
    /* ---------------------------------------------------------------------- */
    /**
     * Nodes coordinate construction and synchronization
     */
    std::multimap< UInt, std::pair<UInt, UInt> > nodes_to_proc;
    /// get the list of nodes to send and send them
    Real * local_nodes = NULL;
    UInt nb_nodes_per_proc[nb_proc];
    UInt * nodes_per_proc[nb_proc];

    comm.broadcast(&(mesh.nb_global_nodes), 1, root);

    /* --------<<<<-NB_NODES + NODES----------------------------------------- */
    for (UInt p = 0; p < nb_proc; ++p) {
      UInt nb_nodes = 0;
      //      UInt * buffer;
      if(p != root) {
	AKANTU_DEBUG_INFO("Receiving number of nodes from proc " << p << " TAG("<< Tag::genTag(p, 0, TAG_NB_NODES) <<")");
	comm.receive(&nb_nodes, 1, p, Tag::genTag(p, 0, TAG_NB_NODES));
	nodes_per_proc[p] = new UInt[nb_nodes];
	nb_nodes_per_proc[p] = nb_nodes;
	AKANTU_DEBUG_INFO("Receiving list of nodes from proc " << p << " TAG("<< Tag::genTag(p, 0, TAG_NODES) <<")");
	comm.receive(nodes_per_proc[p], nb_nodes, p, Tag::genTag(p, 0, TAG_NODES));
      } else {
	nb_nodes = old_nodes->getSize();
	nb_nodes_per_proc[p] = nb_nodes;
	nodes_per_proc[p] = old_nodes->storage();
      }

      /// get the coordinates for the selected nodes
      Real * nodes_to_send = new Real[nb_nodes * spatial_dimension];
      Real * nodes_to_send_tmp = nodes_to_send;
      for (UInt n = 0; n < nb_nodes; ++n) {
	memcpy(nodes_to_send_tmp,
	       nodes->storage() + spatial_dimension * nodes_per_proc[p][n],
	       spatial_dimension * sizeof(Real));
	// nodes_to_proc.insert(std::make_pair(buffer[n], std::make_pair(p, n)));
	nodes_to_send_tmp += spatial_dimension;
      }

      /* -------->>>>-COORDINATES-------------------------------------------- */
      if(p != root) { /// send them for distant processors
	AKANTU_DEBUG_INFO("Sending coordinates to proc " << p << " TAG("<< Tag::genTag(my_rank, 0, TAG_COORDINATES) <<")");
	comm.send(nodes_to_send, nb_nodes * spatial_dimension, p, Tag::genTag(my_rank, 0, TAG_COORDINATES));
	delete [] nodes_to_send;
      } else { /// save them for local processor
	local_nodes = nodes_to_send;
      }
    }


    /// construct the local nodes coordinates
    UInt nb_nodes = old_nodes->getSize();
    nodes->resize(nb_nodes);
    memcpy(nodes->storage(), local_nodes, nb_nodes * spatial_dimension * sizeof(Real));
    delete [] local_nodes;

    Array<Int> * nodes_type_per_proc[nb_proc];
    for (UInt p = 0; p < nb_proc; ++p) {
      nodes_type_per_proc[p] = new Array<Int>(nb_nodes_per_proc[p]);
    }

    communicator.fillNodesType(mesh);

    /* --------<<<<-NODES_TYPE-1--------------------------------------------- */
    for (UInt p = 0; p < nb_proc; ++p) {
      if(p != root) {
	AKANTU_DEBUG_INFO("Receiving first nodes types from proc " << p << " TAG("<< Tag::genTag(my_rank, count, TAG_NODES_TYPE) <<")");
	comm.receive(nodes_type_per_proc[p]->storage(),
		     nb_nodes_per_proc[p], p, Tag::genTag(p, 0, TAG_NODES_TYPE));
      } else {
	nodes_type_per_proc[p]->copy(mesh.getNodesType());
      }
      for (UInt n = 0; n < nb_nodes_per_proc[p]; ++n) {
	if((*nodes_type_per_proc[p])(n) == -2)
	  nodes_to_proc.insert(std::make_pair(nodes_per_proc[p][n], std::make_pair(p, n)));
      }
    }

    std::multimap< UInt, std::pair<UInt, UInt> >::iterator it_node;
    std::pair< std::multimap< UInt, std::pair<UInt, UInt> >::iterator,
	       std::multimap< UInt, std::pair<UInt, UInt> >::iterator > it_range;
    for (UInt i = 0; i < mesh.nb_global_nodes; ++i) {
      it_range = nodes_to_proc.equal_range(i);
      if(it_range.first == nodes_to_proc.end() || it_range.first->first != i) continue;

      UInt node_type = (it_range.first)->second.first;
      for (it_node = it_range.first; it_node != it_range.second; ++it_node) {
	UInt proc = it_node->second.first;
	UInt node = it_node->second.second;
	if(proc != node_type)
	  nodes_type_per_proc[proc]->storage()[node] = node_type;
      }
    }

    /* -------->>>>-NODES_TYPE-2--------------------------------------------- */
    std::vector<CommunicationRequest *> requests;
    for (UInt p = 0; p < nb_proc; ++p) {
      if(p != root) {
	AKANTU_DEBUG_INFO("Sending nodes types to proc " << p << " TAG("<< Tag::genTag(my_rank, 0, TAG_NODES_TYPE) <<")");
	requests.push_back(comm.asyncSend(nodes_type_per_proc[p]->storage(),
					  nb_nodes_per_proc[p], p, Tag::genTag(my_rank, 0, TAG_NODES_TYPE)));
      } else {
	mesh.getNodesTypePointer()->copy(*nodes_type_per_proc[p]);
      }
    }

    comm.waitAll(requests);
    comm.freeCommunicationRequest(requests);
    requests.clear();

    for (UInt p = 0; p < nb_proc; ++p) {
      if(p != root) delete [] nodes_per_proc[p];
      delete nodes_type_per_proc[p];
    }

    /* -------->>>>-NODE GROUPS --------------------------------------------- */
    synchronizeNodeGroupsMaster(communicator, root, mesh);

    /* ---------------------------------------------------------------------- */
    /*  Distant (rank != root)                                                */
    /* ---------------------------------------------------------------------- */
  } else {
    /**
     * connectivity and communications scheme construction on distant processors
     */
    ElementType type = _not_defined;
    UInt count = 0;
    do {
      /* --------<<<<-SIZE--------------------------------------------------- */
      UInt size[5] = { 0 };
      comm.receive(size, 5, root, Tag::genTag(root, count, TAG_SIZES));

      type        = (ElementType) size[0];
      UInt nb_local_element     = size[1];
      UInt nb_ghost_element     = size[2];
      UInt nb_element_to_send   = size[3];
      UInt nb_tags              = size[4];

      if(type != _not_defined) {
	UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
	/* --------<<<<-CONNECTIVITY----------------------------------------- */
	local_connectivity = new UInt[(nb_local_element + nb_ghost_element) *
				      nb_nodes_per_element];
	AKANTU_DEBUG_INFO("Receiving connectivities from proc " << root);
	comm.receive(local_connectivity, nb_nodes_per_element * (nb_local_element +
								  nb_ghost_element),
		     root, Tag::genTag(root, count, TAG_CONNECTIVITY));

	AKANTU_DEBUG_INFO("Renumbering local connectivities");
	MeshUtils::renumberMeshNodes(mesh,
				     local_connectivity,
				     nb_local_element,
				     nb_ghost_element,
				     type,
				     *old_nodes);

	delete [] local_connectivity;

	/* --------<<<<-PARTITIONS--------------------------------------------- */
	local_partitions = new UInt[nb_element_to_send + nb_ghost_element * 2];
	AKANTU_DEBUG_INFO("Receiving partition informations from proc " << root);
	comm.receive(local_partitions,
		     nb_element_to_send + nb_ghost_element * 2,
		     root, Tag::genTag(root, count, TAG_PARTITIONS));

	if(Mesh::getSpatialDimension(type) == mesh.getSpatialDimension()) {
	  AKANTU_DEBUG_INFO("Creating communications scheme");
	  communicator.fillCommunicationScheme(local_partitions,
					       nb_local_element,
					       nb_ghost_element,
					       type);
	}
	delete [] local_partitions;

	/* --------<<<<-TAGS------------------------------------------------- */
	synchronizeTagsRecv(communicator, root, mesh, nb_tags, type,
			    nb_local_element,
			    nb_ghost_element);

	/* --------<<<<-GROUPS----------------------------------------------- */
	synchronizeElementGroups(communicator, root, mesh, type);
      }
      ++count;
    } while(type != _not_defined);

    /**
     * Nodes coordinate construction and synchronization on distant processors
     */
    comm.broadcast(&(mesh.nb_global_nodes), 1, root);

    /* -------->>>>-NB_NODES + NODES----------------------------------------- */
    AKANTU_DEBUG_INFO("Sending list of nodes to proc " << root);
    UInt nb_nodes = old_nodes->getSize();
    comm.send(&nb_nodes, 1, root, Tag::genTag(my_rank, 0, TAG_NB_NODES));
    comm.send(old_nodes->storage(), nb_nodes, root, Tag::genTag(my_rank, 0, TAG_NODES));

    /* --------<<<<-COORDINATES---------------------------------------------- */
    nodes->resize(nb_nodes);
    AKANTU_DEBUG_INFO("Receiving coordinates from proc " << root);
    comm.receive(nodes->storage(), nb_nodes * spatial_dimension, root, Tag::genTag(root, 0, TAG_COORDINATES));

    communicator.fillNodesType(mesh);
    /* -------->>>>-NODES_TYPE-1--------------------------------------------- */
    Int * nodes_types = mesh.getNodesTypePointer()->storage();
    AKANTU_DEBUG_INFO("Sending first nodes types to proc " << root);
    comm.send(nodes_types, nb_nodes,
	      root, Tag::genTag(my_rank, 0, TAG_NODES_TYPE));

    /* --------<<<<-NODES_TYPE-2--------------------------------------------- */
    AKANTU_DEBUG_INFO("Receiving nodes types from proc " << root);
    comm.receive(nodes_types, nb_nodes,
		 root, Tag::genTag(root, 0, TAG_NODES_TYPE));

    /* --------<<<<-NODE GROUPS --------------------------------------------- */
    synchronizeNodeGroupsSlaves(communicator, root, mesh);
  }

  MeshUtils::fillElementToSubElementsData(mesh);

  mesh.is_distributed = true;

  AKANTU_DEBUG_OUT();
  return &communicator;
}

/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::fillTagBuffer(const MeshData & mesh_data,
					    DynamicCommunicationBuffer * buffers,
					    const std::string & tag_name,
					    const ElementType & el_type,
					    const Array<UInt> & partition_num,
					    const CSR<UInt> & ghost_partition) {
  #define AKANTU_DISTRIBUTED_SYNHRONIZER_TAG_DATA(r, extra_param, elem)	\
    case BOOST_PP_TUPLE_ELEM(2, 0, elem) : { \
      fillTagBufferTemplated<BOOST_PP_TUPLE_ELEM(2, 1, elem)>(mesh_data, buffers, tag_name, el_type, partition_num, ghost_partition); \
      break; \
    } \

  MeshDataTypeCode data_type_code = mesh_data.getTypeCode(tag_name);
  switch(data_type_code) {
    BOOST_PP_SEQ_FOR_EACH(AKANTU_DISTRIBUTED_SYNHRONIZER_TAG_DATA, , AKANTU_MESH_DATA_TYPES)
  default : AKANTU_DEBUG_ERROR("Could not obtain the type of tag" << tag_name << "!"); break;
  }
  #undef AKANTU_DISTRIBUTED_SYNHRONIZER_TAG_DATA
}

/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::fillNodesType(Mesh & mesh) {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = mesh.getNbNodes();
  Int * nodes_type = mesh.getNodesTypePointer()->storage();

  UInt * nodes_set = new UInt[nb_nodes];
  std::fill_n(nodes_set, nb_nodes, 0);

  const UInt NORMAL_SET = 1;
  const UInt GHOST_SET  = 2;

  bool * already_seen = new bool[nb_nodes];

  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;
    UInt set = NORMAL_SET;
    if (gt == _ghost) set = GHOST_SET;

    std::fill_n(already_seen, nb_nodes, false);
    Mesh::type_iterator it  = mesh.firstType(_all_dimensions, gt, _ek_not_defined);
    Mesh::type_iterator end = mesh.lastType(_all_dimensions, gt, _ek_not_defined);
    for(; it != end; ++it) {
      ElementType type = *it;

      UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
      UInt nb_element = mesh.getNbElement(type, gt);
      Array<UInt>::iterator< Vector<UInt> > conn_it = mesh.getConnectivity(type, gt).begin(nb_nodes_per_element);

      for (UInt e = 0; e < nb_element; ++e, ++conn_it) {
	Vector<UInt> & conn = *conn_it;
	for (UInt n = 0; n < nb_nodes_per_element; ++n) {
	  AKANTU_DEBUG_ASSERT(conn(n) < nb_nodes, "Node " << conn(n)
			      << " bigger than number of nodes " << nb_nodes);
	  if(!already_seen[conn(n)]) {
	    nodes_set[conn(n)] += set;
	    already_seen[conn(n)] = true;
	  }
	}
      }
    }
  }

  delete [] already_seen;

  for (UInt i = 0; i < nb_nodes; ++i) {
    if(nodes_set[i] == NORMAL_SET) nodes_type[i] = -1;
    else if(nodes_set[i] == GHOST_SET) nodes_type[i] = -3;
    else if(nodes_set[i] == (GHOST_SET + NORMAL_SET)) nodes_type[i] = -2;
  }

  delete [] nodes_set;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::fillCommunicationScheme(const UInt * partition,
						      UInt nb_local_element,
						      UInt nb_ghost_element,
						      ElementType type) {
  AKANTU_DEBUG_IN();

  Element element;
  element.type = type;
  element.kind = Mesh::getKind(type);

  const UInt * part = partition;

  part = partition;
  for (UInt lel = 0; lel < nb_local_element; ++lel) {
    UInt nb_send = *part; part++;
    element.element = lel;
    element.ghost_type = _not_ghost;
    for (UInt p = 0; p < nb_send; ++p) {
      UInt proc = *part; part++;

      AKANTU_DEBUG(dblAccessory, "Must send : " << element << " to proc " << proc);
      (send_element[proc]).push_back(element);
    }
  }

  for (UInt gel = 0; gel < nb_ghost_element; ++gel) {
    UInt proc = *part; part++;
    element.element = gel;
    element.ghost_type = _ghost;
    AKANTU_DEBUG(dblAccessory, "Must recv : " << element << " from proc " << proc);
    recv_element[proc].push_back(element);
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::asynchronousSynchronize(DataAccessor & data_accessor,
						      SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  if (communications.find(tag) == communications.end())
    computeBufferSize(data_accessor, tag);

  Communication & communication = communications[tag];

  AKANTU_DEBUG_ASSERT(communication.send_requests.size() == 0,
		      "There must be some pending sending communications. Tag is " << tag);

  std::map<SynchronizationTag, UInt>::iterator t_it = tag_counter.find(tag);
  UInt counter = 0;
  if(t_it == tag_counter.end()) {
    tag_counter[tag] = 0;
  } else {
    counter = ++(t_it->second);
  }

  for (UInt p = 0; p < nb_proc; ++p) {
    UInt ssize = communication.size_to_send[p];
    if(p == rank || ssize == 0) continue;

    CommunicationBuffer & buffer = communication.send_buffer[p];
    buffer.resize(ssize);
#ifndef AKANTU_NDEBUG
    UInt nb_elements   =  send_element[p].getSize();
    AKANTU_DEBUG_INFO("Packing data for proc " << p
		      << " (" << ssize << "/" << nb_elements
		      <<" data to send/elements)");

    /// pack barycenters in debug mode
    Array<Element>::const_iterator<Element> bit  = send_element[p].begin();
    Array<Element>::const_iterator<Element> bend = send_element[p].end();
    for (; bit != bend; ++bit) {
      const Element & element = *bit;
      Vector<Real> barycenter(mesh.getSpatialDimension());
      mesh.getBarycenter(element.element, element.type, barycenter.storage(), element.ghost_type);
      buffer << barycenter;
    }
#endif

    data_accessor.packElementData(buffer, send_element[p], tag);

    AKANTU_DEBUG_ASSERT(buffer.getPackedSize() == ssize,
			"a problem have been introduced with "
			<< "false sent sizes declaration "
			<< buffer.getPackedSize() << " != " << ssize);
    AKANTU_DEBUG_INFO("Posting send to proc " << p
		      << " (tag: " << tag << " - " << ssize << " data to send)"
		      << " [" << Tag::genTag(rank, counter, tag) << "]");
    communication.send_requests.push_back(static_communicator->asyncSend(buffer.storage(),
									 ssize,
									 p,
									 Tag::genTag(rank, counter, tag)));
  }

  AKANTU_DEBUG_ASSERT(communication.recv_requests.size() == 0,
		      "There must be some pending receive communications");

  for (UInt p = 0; p < nb_proc; ++p) {
    UInt rsize = communication.size_to_receive[p];
    if(p == rank || rsize == 0) continue;
    CommunicationBuffer & buffer = communication.recv_buffer[p];
    buffer.resize(rsize);

    AKANTU_DEBUG_INFO("Posting receive from proc " << p
		      << " (tag: " << tag << " - " << rsize << " data to receive) "
		      << " [" << Tag::genTag(p, counter, tag) << "]");
    communication.recv_requests.push_back(static_communicator->asyncReceive(buffer.storage(),
									    rsize,
									    p,
									    Tag::genTag(p, counter, tag)));
  }


#if defined(AKANTU_DEBUG_TOOLS) && defined(AKANTU_CORE_CXX11)
  static std::set<SynchronizationTag> tags;
  if(tags.find(tag) == tags.end()) {
    debug::element_manager.print(debug::_dm_synch,
				 [&send_element, rank, nb_proc, tag, id](const Element & el)->std::string {
				   std::stringstream out;
				   UInt elp = 0;
				   for (UInt p = 0; p < nb_proc; ++p) {
				     UInt pos = send_element[p].find(el);
				     if(pos != UInt(-1)) {
				       if(elp > 0) out << std::endl;
				       out << id << " send (" << pos << "/" << send_element[p].getSize() << ") to proc " << p << " tag:" << tag;
				       ++elp;
				     }
				   }
				   return out.str();
				 });

    debug::element_manager.print(debug::_dm_synch,
				 [&recv_element, rank, nb_proc, tag, id](const Element & el)->std::string {
				   std::stringstream out;
				   UInt elp = 0;
				   for (UInt p = 0; p < nb_proc; ++p) {
				     if(p == rank) continue;
				     UInt pos = recv_element[p].find(el);
				     if(pos != UInt(-1)) {
				       if(elp > 0) out << std::endl;
				       out << id << " recv (" << pos << "/" << recv_element[p].getSize() << ") from proc " << p << " tag:" << tag;
				       ++elp;
				     }
				   }
				   return out.str();
				 });
    tags.insert(tag);
  }
#endif


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::waitEndSynchronize(DataAccessor & data_accessor,
						 SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(communications.find(tag) != communications.end(), "No communication with the tag \""
		      << tag <<"\" started");

  Communication & communication = communications[tag];

  std::vector<CommunicationRequest *> req_not_finished;
  std::vector<CommunicationRequest *> * req_not_finished_tmp = &req_not_finished;
  std::vector<CommunicationRequest *> * recv_requests_tmp = &(communication.recv_requests);

  //  static_communicator->waitAll(recv_requests);
  while(!recv_requests_tmp->empty()) {
    for (std::vector<CommunicationRequest *>::iterator req_it = recv_requests_tmp->begin();
	 req_it != recv_requests_tmp->end() ; ++req_it) {
      CommunicationRequest * req = *req_it;

      if(static_communicator->testRequest(req)) {
	UInt proc = req->getSource();
	AKANTU_DEBUG_INFO("Unpacking data coming from proc " << proc);
	CommunicationBuffer & buffer = communication.recv_buffer[proc];

#ifndef AKANTU_NDEBUG
	Array<Element>::const_iterator<Element> bit  = recv_element[proc].begin();
	Array<Element>::const_iterator<Element> bend = recv_element[proc].end();

	UInt spatial_dimension = mesh.getSpatialDimension();

	for (; bit != bend; ++bit) {
	  const Element & element = *bit;

	  Vector<Real> barycenter_loc(spatial_dimension);
	  mesh.getBarycenter(element.element,
			     element.type,
			     barycenter_loc.storage(),
			     element.ghost_type);
	  Vector<Real> barycenter(spatial_dimension);
	  buffer >> barycenter;
	  Real tolerance = Math::getTolerance();
	  Real bary_norm = barycenter.norm();
	  for (UInt i = 0; i < spatial_dimension; ++i) {
	    if((std::abs((barycenter(i) - barycenter_loc(i))/bary_norm) <= tolerance) ||
	       (std::abs(barycenter_loc(i)) <= 0 && std::abs(barycenter(i)) <= tolerance)) continue;
	    AKANTU_DEBUG_ERROR("Unpacking an unknown value for the element: "
			       << element
			       << "(barycenter[" << i << "] = " << barycenter_loc(i)
			       << " and buffer[" << i << "] = " << barycenter(i) << ") ["
			       << std::abs((barycenter(i) - barycenter_loc(i))/barycenter_loc(i))
			       << "] - tag: " << tag);
	  }
	}
#endif

	data_accessor.unpackElementData(buffer, recv_element[proc], tag);
	buffer.resize(0);

	AKANTU_DEBUG_ASSERT(buffer.getLeftToUnpack() == 0,
			    "all data have not been unpacked: "
			    << buffer.getLeftToUnpack() << " bytes left");
	static_communicator->freeCommunicationRequest(req);
      } else {
	req_not_finished_tmp->push_back(req);
      }
    }

    std::vector<CommunicationRequest *> * swap = req_not_finished_tmp;
    req_not_finished_tmp = recv_requests_tmp;
    recv_requests_tmp = swap;

    req_not_finished_tmp->clear();
  }


  AKANTU_DEBUG_INFO("Waiting that every send requests are received");
  static_communicator->waitAll(communication.send_requests);
  for (std::vector<CommunicationRequest *>::iterator req_it = communication.send_requests.begin();
       req_it != communication.send_requests.end() ; ++req_it) {
    CommunicationRequest & req = *(*req_it);

    if(static_communicator->testRequest(&req)) {
      UInt proc = req.getDestination();
      CommunicationBuffer & buffer = communication.send_buffer[proc];
      buffer.resize(0);
      static_communicator->freeCommunicationRequest(&req);
    }
  }
  communication.send_requests.clear();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::computeBufferSize(DataAccessor & data_accessor,
						SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  communications[tag].resize(nb_proc);

  for (UInt p = 0; p < nb_proc; ++p) {
    UInt ssend    = 0;
    UInt sreceive = 0;
    if(p != rank) {
      if(send_element[p].getSize() != 0) {
#ifndef AKANTU_NDEBUG
	ssend += send_element[p].getSize() * mesh.getSpatialDimension() * sizeof(Real);
#endif
	ssend += data_accessor.getNbDataForElements(send_element[p], tag);
	AKANTU_DEBUG_INFO("I have " << ssend << "(" << ssend / 1024.
			  << "kB - "<< send_element[p].getSize() <<" element(s)) data to send to " << p << " for tag "
			  << tag);
      }

      if(recv_element[p].getSize() != 0) {
#ifndef AKANTU_NDEBUG
	sreceive += recv_element[p].getSize() * mesh.getSpatialDimension() * sizeof(Real);
#endif
	sreceive += data_accessor.getNbDataForElements(recv_element[p], tag);
	AKANTU_DEBUG_INFO("I have " << sreceive << "(" << sreceive / 1024.
			  << "kB - "<< recv_element[p].getSize() <<" element(s)) data to receive for tag "
			  << tag);
      }
    }

    communications[tag].size_to_send   [p] = ssend;
    communications[tag].size_to_receive[p] = sreceive;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  Int prank = StaticCommunicator::getStaticCommunicator().whoAmI();
  Int psize = StaticCommunicator::getStaticCommunicator().getNbProc();
  stream << "[" << prank << "/" << psize << "]" << space << "DistributedSynchronizer [" << std::endl;
  for (UInt p = 0; p < nb_proc; ++p) {
    if (p == UInt(prank)) continue;
    stream << "[" << prank << "/" << psize << "]" << space
	   << " + Communication to proc " << p << " [" << std::endl;
    if(AKANTU_DEBUG_TEST(dblDump)) {
      stream << "[" << prank << "/" << psize << "]" << space
	     << "    - Element to send to proc " << p << " [" << std::endl;

      Array<Element>::iterator<Element> it_el  = send_element[p].begin();
      Array<Element>::iterator<Element> end_el = send_element[p].end();
      for(;it_el != end_el; ++it_el)
	stream << "[" << prank << "/" << psize << "]" << space << "       " << *it_el << std::endl;
      stream << "[" << prank << "/" << psize << "]" << space << "   ]" << std::endl;

      stream << "[" << prank << "/" << psize << "]" << space
	     << "    - Element to recv from proc " << p << " [" << std::endl;

      it_el  = recv_element[p].begin();
      end_el = recv_element[p].end();
      for(;it_el != end_el; ++it_el)
	stream << "[" << prank << "/" << psize << "]"
	       << space << "       " << *it_el << std::endl;

      stream << "[" << prank << "/" << psize << "]" << space << "   ]" << std::endl;
    }

    std::map< SynchronizationTag, Communication>::const_iterator it = communications.begin();
    std::map< SynchronizationTag, Communication>::const_iterator end = communications.end();
    for (; it != end; ++it) {
      const SynchronizationTag & tag = it->first;
      const Communication & communication = it->second;
      UInt ssend    = communication.size_to_send[p];
      UInt sreceive = communication.size_to_receive[p];
      stream << "[" << prank << "/" << psize << "]" << space << "     - Tag " << tag << " -> " << ssend << "byte(s) -- <- " << sreceive << "byte(s)" << std::endl;
    }
  }
}

/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::onElementsRemoved(const Array<Element> & element_to_remove,
						const ElementTypeMapArray<UInt> & new_numbering,
						__attribute__((unused)) const RemovedElementsEvent & event) {
  AKANTU_DEBUG_IN();

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  UInt psize = comm.getNbProc();
  UInt prank = comm.whoAmI();

  std::vector<CommunicationRequest *> isend_requests;
  Array<UInt> * list_of_el = new Array<UInt>[nb_proc];
  // Handling ghost elements
  for (UInt p = 0; p < psize; ++p) {
    if (p == prank) continue;

    Array<Element> & recv = recv_element[p];
    if(recv.getSize() == 0) continue;

    Array<Element>::iterator<Element> recv_begin = recv.begin();
    Array<Element>::iterator<Element> recv_end   = recv.end();

    Array<Element>::const_iterator<Element> er_it  = element_to_remove.begin();
    Array<Element>::const_iterator<Element> er_end = element_to_remove.end();

    Array<UInt> & list = list_of_el[p];
    for (UInt i = 0; recv_begin != recv_end; ++i, ++recv_begin) {
      const Element & el = *recv_begin;
      Array<Element>::const_iterator<Element> pos = std::find(er_it, er_end, el);
      if(pos == er_end) {
	list.push_back(i);
      }
    }

    if(list.getSize() == recv.getSize())
      list.push_back(UInt(0));
    else list.push_back(UInt(-1));

    AKANTU_DEBUG_INFO("Sending a message of size " << list.getSize() << " to proc " << p << " TAG(" << Tag::genTag(prank, 0, 0) << ")");
    isend_requests.push_back(comm.asyncSend(list.storage(), list.getSize(),
					    p, Tag::genTag(prank, 0, 0)));

    list.erase(list.getSize() - 1);
    if(list.getSize() == recv.getSize()) continue;

    Array<Element> new_recv;
    for (UInt nr = 0; nr < list.getSize(); ++nr) {
      Element & el = recv(list(nr));
      el.element = new_numbering(el.type, el.ghost_type)(el.element);
      new_recv.push_back(el);
    }

    AKANTU_DEBUG_INFO("I had " << recv.getSize() << " elements to recv from proc " << p << " and "
		      << list.getSize() << " elements to keep. I have "
		      << new_recv.getSize() << " elements left.");
    recv.copy(new_recv);
  }

  for (UInt p = 0; p < psize; ++p) {
    if (p == prank) continue;
    Array<Element> & send = send_element[p];

    if(send.getSize() == 0) continue;

    CommunicationStatus status;
    AKANTU_DEBUG_INFO("Getting number of elements of proc " << p << " not needed anymore TAG("<< Tag::genTag(p, 0, 0) <<")");
    comm.probe<UInt>(p, Tag::genTag(p, 0, 0), status);
    Array<UInt> list(status.getSize());

    AKANTU_DEBUG_INFO("Receiving list of elements (" << status.getSize() - 1
		      << " elements) no longer needed by proc " << p
		      << " TAG("<< Tag::genTag(p, 0, 0) <<")");
    comm.receive(list.storage(), list.getSize(),
		 p, Tag::genTag(p, 0, 0));

    if(list.getSize() == 1 && list(0) == 0) continue;

    list.erase(list.getSize() - 1);

    Array<Element> new_send;
    for (UInt ns = 0; ns < list.getSize(); ++ns) {
      new_send.push_back(send(list(ns)));
    }

    AKANTU_DEBUG_INFO("I had " << send.getSize() << " elements to send to proc " << p << " and "
		      << list.getSize() << " elements to keep. I have "
		      << new_send.getSize() << " elements left.");
    send.copy(new_send);
  }

  comm.waitAll(isend_requests);
  comm.freeCommunicationRequest(isend_requests);

  delete [] list_of_el;
  AKANTU_DEBUG_OUT();
}


// void DistributedSynchronizer::checkCommunicationScheme() {
//   for (UInt p = 0; p < psize; ++p) {
//     if (p == prank) continue;
//     for(UInt e(0), e < recv_element.getSize())
//   }
// }

/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::buildPrankToElement() {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();

  mesh.initElementTypeMapArray(prank_to_element,
			      1,
			      spatial_dimension,
			      false,
			      _ek_not_defined,
			      true);

  Mesh::type_iterator it  = mesh.firstType(spatial_dimension,
					   _not_ghost,
					   _ek_not_defined);

  Mesh::type_iterator end = mesh.lastType(spatial_dimension,
					  _not_ghost,
					  _ek_not_defined);

  /// assign prank to all not ghost elements
  for (; it != end; ++it) {
    UInt nb_element = mesh.getNbElement(*it);
    Array<UInt> & prank_to_el = prank_to_element(*it);
    for (UInt el = 0; el < nb_element; ++el) {
      prank_to_el(el) = rank;
    }
  }

  /// assign prank to all ghost elements
  for (UInt p = 0; p < nb_proc; ++p) {
    UInt nb_ghost_element = recv_element[p].getSize();

    for (UInt el = 0; el < nb_ghost_element; ++el) {
      UInt element = recv_element[p](el).element;
      ElementType type = recv_element[p](el).type;
      GhostType ghost_type = recv_element[p](el).ghost_type;

      Array<UInt> & prank_to_el = prank_to_element(type, ghost_type);
      prank_to_el(element) = p;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::filterElementsByKind(DistributedSynchronizer * new_synchronizer,
						   ElementKind kind) {
  AKANTU_DEBUG_IN();

  Array<Element> * newsy_send_element = new_synchronizer->send_element;
  Array<Element> * newsy_recv_element = new_synchronizer->recv_element;

  Array<Element> * new_send_element = new Array<Element>[nb_proc];
  Array<Element> * new_recv_element = new Array<Element>[nb_proc];

  for (UInt p = 0; p < nb_proc; ++p) {

    /// send element copying part
    new_send_element[p].resize(0);

    for (UInt el = 0; el < send_element[p].getSize(); ++el) {
      Element & element = send_element[p](el);

      if (element.kind == kind)
	newsy_send_element[p].push_back(element);
      else
	new_send_element[p].push_back(element);
    }

    /// recv element copying part
    new_recv_element[p].resize(0);

    for (UInt el = 0; el < recv_element[p].getSize(); ++el) {
      Element & element = recv_element[p](el);

      if (element.kind == kind)
	newsy_recv_element[p].push_back(element);
      else
	new_recv_element[p].push_back(element);
    }
  }

  /// deleting and reassigning old pointers
  delete [] send_element;
  delete [] recv_element;

  send_element = new_send_element;
  recv_element = new_recv_element;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::reset() {
  AKANTU_DEBUG_IN();

  for (UInt p = 0; p < nb_proc; ++p) {
    send_element[p].resize(0);
    recv_element[p].resize(0);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::synchronizeTagsSend(DistributedSynchronizer & communicator,
						  UInt root,
						  Mesh & mesh,
						  UInt nb_tags,
						  const ElementType & type,
						  const Array<UInt> & partition_num,
						  const CSR<UInt> & ghost_partition,
						  UInt nb_local_element,
						  UInt nb_ghost_element) {
  AKANTU_DEBUG_IN();

  static UInt count = 0;

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  UInt nb_proc = comm.getNbProc();
  UInt my_rank = comm.whoAmI();

  if(nb_tags == 0) {
    AKANTU_DEBUG_OUT();
    return;
  }

  UInt mesh_data_sizes_buffer_length;
  MeshData & mesh_data = mesh.getMeshData();

  /// tag info
  std::vector<std::string> tag_names;
  mesh.getMeshData().getTagNames(tag_names, type);
  // Make sure the tags are sorted (or at least not in random order),
  // because they come from a map !!
  std::sort(tag_names.begin(), tag_names.end());


  // Sending information about the tags in mesh_data: name, data type and
  // number of components of the underlying array associated to the current type
  DynamicCommunicationBuffer mesh_data_sizes_buffer;
  std::vector<std::string>::const_iterator names_it  = tag_names.begin();
  std::vector<std::string>::const_iterator names_end = tag_names.end();
  for(;names_it != names_end; ++names_it) {
    mesh_data_sizes_buffer << *names_it;
    mesh_data_sizes_buffer << mesh_data.getTypeCode(*names_it);
    mesh_data_sizes_buffer << mesh_data.getNbComponent(*names_it, type);
  }

  mesh_data_sizes_buffer_length = mesh_data_sizes_buffer.getSize();
  AKANTU_DEBUG_INFO("Broadcasting the size of the information about the mesh data tags: (" << mesh_data_sizes_buffer_length << ")." );
  comm.broadcast(&mesh_data_sizes_buffer_length, 1, root);
  AKANTU_DEBUG_INFO("Broadcasting the information about the mesh data tags, addr " << (void*)mesh_data_sizes_buffer.storage());

  if(mesh_data_sizes_buffer_length !=0)
    comm.broadcast(mesh_data_sizes_buffer.storage(), mesh_data_sizes_buffer.getSize(), root);


  if(mesh_data_sizes_buffer_length !=0) {
    //Sending the actual data to each processor
    DynamicCommunicationBuffer buffers[nb_proc];
    std::vector<std::string>::const_iterator names_it  = tag_names.begin();
    std::vector<std::string>::const_iterator names_end = tag_names.end();

    // Loop over each tag for the current type
    for(;names_it != names_end; ++names_it) {
      // Type code of the current tag (i.e. the tag named *names_it)
      communicator.fillTagBuffer(mesh_data,
				 buffers,
				 *names_it,
				 type,
				 partition_num,
				 ghost_partition);
    }

    std::vector<CommunicationRequest *> requests;
    for (UInt p = 0; p < nb_proc; ++p) {
      if(p != root) {
	AKANTU_DEBUG_INFO("Sending " << buffers[p].getSize() << " bytes of mesh data to proc " << p << " TAG("<< Tag::genTag(my_rank, count, TAG_MESH_DATA) <<")");

	requests.push_back(comm.asyncSend(buffers[p].storage(),
					  buffers[p].getSize(), p, Tag::genTag(my_rank, count, TAG_MESH_DATA)));
      }
    }

    names_it  = tag_names.begin();
    // Loop over each tag for the current type
    for(;names_it != names_end; ++names_it) {
      // Reinitializing the mesh data on the master
      communicator.populateMeshData(mesh_data,
				    buffers[root],
				    *names_it,
				    type,
				    mesh_data.getTypeCode(*names_it),
				    mesh_data.getNbComponent(*names_it, type),
				    nb_local_element,
				    nb_ghost_element);
    }


    comm.waitAll(requests);
    comm.freeCommunicationRequest(requests);
    requests.clear();
  }

  ++count;

  AKANTU_DEBUG_OUT();
}



/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::synchronizeTagsRecv(DistributedSynchronizer & communicator,
						  UInt root,
						  Mesh & mesh,
						  UInt nb_tags,
						  const ElementType & type,
						  UInt nb_local_element,
						  UInt nb_ghost_element) {
  AKANTU_DEBUG_IN();

  static UInt count = 0;

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();

  if(nb_tags == 0) {
    AKANTU_DEBUG_OUT();
    return;
  }

  /* --------<<<<-TAGS------------------------------------------------- */
  UInt mesh_data_sizes_buffer_length = 0;
  CommunicationBuffer mesh_data_sizes_buffer;
  MeshData & mesh_data = mesh.getMeshData();

  AKANTU_DEBUG_INFO("Receiving the size of the information about the mesh data tags.");
  comm.broadcast(&mesh_data_sizes_buffer_length, 1, root);

  if(mesh_data_sizes_buffer_length != 0) {
    mesh_data_sizes_buffer.resize(mesh_data_sizes_buffer_length);
    AKANTU_DEBUG_INFO("Receiving the information about the mesh data tags, addr " << (void*)mesh_data_sizes_buffer.storage());
    comm.broadcast(mesh_data_sizes_buffer.storage(), mesh_data_sizes_buffer_length, root);
    AKANTU_DEBUG_INFO("Size of the information about the mesh data: " << mesh_data_sizes_buffer_length);

    std::vector<std::string> tag_names;
    std::vector<MeshDataTypeCode> tag_type_codes;
    std::vector<UInt> tag_nb_component;
    tag_names.resize(nb_tags);
    tag_type_codes.resize(nb_tags);
    tag_nb_component.resize(nb_tags);
    CommunicationBuffer mesh_data_buffer;
    UInt type_code_int;
    for(UInt i(0); i < nb_tags; ++i) {
      mesh_data_sizes_buffer >> tag_names[i];
      mesh_data_sizes_buffer >> type_code_int;
      tag_type_codes[i] = static_cast<MeshDataTypeCode>(type_code_int);
      mesh_data_sizes_buffer >> tag_nb_component[i];
    }

    std::vector<std::string>::const_iterator names_it  = tag_names.begin();
    std::vector<std::string>::const_iterator names_end = tag_names.end();

    CommunicationStatus mesh_data_comm_status;
    AKANTU_DEBUG_INFO("Checking size of data to receive for mesh data TAG("
		      << Tag::genTag(root, count, TAG_MESH_DATA) << ")");
    comm.probe<char>(root, Tag::genTag(root, count, TAG_MESH_DATA), mesh_data_comm_status);
    UInt mesh_data_buffer_size(mesh_data_comm_status.getSize());
    AKANTU_DEBUG_INFO("Receiving " << mesh_data_buffer_size
		      << " bytes of mesh data TAG("
		      << Tag::genTag(root, count, TAG_MESH_DATA) << ")");
    mesh_data_buffer.resize(mesh_data_buffer_size);
    comm.receive(mesh_data_buffer.storage(),
		 mesh_data_buffer_size,
		 root,
		 Tag::genTag(root, count, TAG_MESH_DATA));

    // Loop over each tag for the current type
    UInt k(0);
    for(; names_it != names_end; ++names_it, ++k) {
      communicator.populateMeshData(mesh_data, mesh_data_buffer,
				    *names_it, type,
				    tag_type_codes[k],
				    tag_nb_component[k],
				    nb_local_element, nb_ghost_element);
    }
  }

  ++count;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<class CommunicationBuffer>
void DistributedSynchronizer::fillElementGroupsFromBuffer(DistributedSynchronizer & communicator,
							  Mesh & mesh,
							  const ElementType & type,
							  CommunicationBuffer & buffer) {
  AKANTU_DEBUG_IN();

  Element el;
  el.type = type;

  for (ghost_type_t::iterator gt = ghost_type_t::begin();  gt != ghost_type_t::end(); ++gt) {
    UInt nb_element = mesh.getNbElement(type, *gt);
    el.ghost_type = *gt;

    for (UInt e = 0; e < nb_element; ++e) {
      el.element = e;

      std::vector<std::string> element_to_group;
      buffer >> element_to_group;

      AKANTU_DEBUG_ASSERT(e < mesh.getNbElement(type, *gt), "The mesh does not have the element " << e);

      std::vector<std::string>::iterator it  = element_to_group.begin();
      std::vector<std::string>::iterator end = element_to_group.end();
      for (; it != end; ++it) {
	mesh.getElementGroup(*it).add(el, false, false);
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::synchronizeElementGroups(DistributedSynchronizer & communicator,
						       UInt root,
						       Mesh & mesh,
						       const ElementType & type,
						       const Array<UInt> & partition_num,
						       const CSR<UInt> & ghost_partition,
						       UInt nb_element) {
  AKANTU_DEBUG_IN();

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();

  UInt nb_proc = comm.getNbProc();
  UInt my_rank = comm.whoAmI();

  DynamicCommunicationBuffer buffers[nb_proc];

  typedef std::vector< std::vector<std::string> > ElementToGroup;
  ElementToGroup element_to_group;
  element_to_group.resize(nb_element);

  GroupManager::const_element_group_iterator egi = mesh.element_group_begin();
  GroupManager::const_element_group_iterator ege = mesh.element_group_end();
  for (; egi != ege; ++egi) {
    ElementGroup & eg = *(egi->second);

    std::string name = egi->first;

    ElementGroup::const_element_iterator eit  = eg.element_begin(type, _not_ghost);
    ElementGroup::const_element_iterator eend = eg.element_end(type, _not_ghost);
    for (; eit != eend; ++eit) {
      element_to_group[*eit].push_back(name);
    }

    eit  = eg.element_begin(type, _not_ghost);
    if(eit != eend)
      const_cast<Array<UInt> &>(eg.getElements(type)).empty();
  }

  /// preparing the buffers
  const UInt * part = partition_num.storage();

  /// copying the data, element by element
  ElementToGroup::const_iterator data_it  = element_to_group.begin();
  ElementToGroup::const_iterator data_end = element_to_group.end();
  for (; data_it != data_end; ++part, ++data_it) {
    buffers[*part] << *data_it;
  }

  data_it = element_to_group.begin();
  /// copying the data for the ghost element
  for (UInt el(0); data_it != data_end; ++data_it, ++el) {
    CSR<UInt>::const_iterator it  = ghost_partition.begin(el);
    CSR<UInt>::const_iterator end = ghost_partition.end(el);
    for (;it != end; ++it) {
      UInt proc = *it;
      buffers[proc] << *data_it;
    }
  }


  std::vector<CommunicationRequest *> requests;
  for (UInt p = 0; p < nb_proc; ++p) {
    if(p == my_rank) continue;
    AKANTU_DEBUG_INFO("Sending element groups to proc " << p
		      << " TAG("<< Tag::genTag(my_rank, p, TAG_ELEMENT_GROUP) <<")");
    requests.push_back(comm.asyncSend(buffers[p].storage(),
				      buffers[p].getSize(),
				      p,
				      Tag::genTag(my_rank, p, TAG_ELEMENT_GROUP)));
  }

  fillElementGroupsFromBuffer(communicator, mesh, type, buffers[my_rank]);

  comm.waitAll(requests);
  comm.freeCommunicationRequest(requests);
  requests.clear();

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::synchronizeElementGroups(DistributedSynchronizer & communicator,
						       UInt root,
						       Mesh & mesh,
						       const ElementType & type) {
  AKANTU_DEBUG_IN();

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  UInt my_rank = comm.whoAmI();


  AKANTU_DEBUG_INFO("Receiving element groups from proc " << root
		    << " TAG("<< Tag::genTag(root, my_rank, TAG_ELEMENT_GROUP) <<")");

  CommunicationStatus status;
  comm.probe<char>(root, Tag::genTag(root, my_rank, TAG_ELEMENT_GROUP), status);

  CommunicationBuffer buffer(status.getSize());
  comm.receive(buffer.storage(), buffer.getSize(), root, Tag::genTag(root, my_rank, TAG_ELEMENT_GROUP));

  fillElementGroupsFromBuffer(communicator, mesh, type, buffer);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<class CommunicationBuffer>
void DistributedSynchronizer::fillNodeGroupsFromBuffer(DistributedSynchronizer & communicator,
						       Mesh & mesh,
						       CommunicationBuffer & buffer) {
  AKANTU_DEBUG_IN();

  std::vector< std::vector<std::string> > node_to_group;

  buffer >> node_to_group;

  AKANTU_DEBUG_ASSERT(node_to_group.size() == mesh.getNbGlobalNodes(),
		      "Not the good amount of nodes where transmitted");

  const Array<UInt> & global_nodes = mesh.getGlobalNodesIds();

  Array<UInt>::const_scalar_iterator nbegin = global_nodes.begin();
  Array<UInt>::const_scalar_iterator nit    = global_nodes.begin();
  Array<UInt>::const_scalar_iterator nend   = global_nodes.end();

  for (; nit != nend; ++nit) {
    std::vector<std::string>::iterator it    = node_to_group[*nit].begin();
    std::vector<std::string>::iterator end   = node_to_group[*nit].end();

    for (; it != end; ++it) {
      mesh.getNodeGroup(*it).add(nit - nbegin, false);
    }
  }

  GroupManager::const_node_group_iterator ngi = mesh.node_group_begin();
  GroupManager::const_node_group_iterator nge = mesh.node_group_end();
  for (; ngi != nge; ++ngi) {
    NodeGroup & ng = *(ngi->second);
    ng.optimize();
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::synchronizeNodeGroupsMaster(DistributedSynchronizer & communicator,
							  UInt root,
							  Mesh & mesh) {
  AKANTU_DEBUG_IN();

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();

  UInt nb_proc = comm.getNbProc();
  UInt my_rank = comm.whoAmI();

  UInt nb_total_nodes = mesh.getNbGlobalNodes();

  DynamicCommunicationBuffer buffer;

  typedef std::vector< std::vector<std::string> > NodeToGroup;
  NodeToGroup node_to_group;
  node_to_group.resize(nb_total_nodes);

  GroupManager::const_node_group_iterator ngi = mesh.node_group_begin();
  GroupManager::const_node_group_iterator nge = mesh.node_group_end();
  for (; ngi != nge; ++ngi) {
    NodeGroup & ng = *(ngi->second);

    std::string name = ngi->first;

    NodeGroup::const_node_iterator nit  = ng.begin();
    NodeGroup::const_node_iterator nend = ng.end();
    for (; nit != nend; ++nit) {
      node_to_group[*nit].push_back(name);
    }

    nit  = ng.begin();
    if(nit != nend)
      ng.empty();
  }

  buffer << node_to_group;

  std::vector<CommunicationRequest *> requests;
  for (UInt p = 0; p < nb_proc; ++p) {
    if(p == my_rank) continue;
    AKANTU_DEBUG_INFO("Sending node groups to proc " << p
		      << " TAG("<< Tag::genTag(my_rank, p, TAG_NODE_GROUP) <<")");
    requests.push_back(comm.asyncSend(buffer.storage(),
				      buffer.getSize(),
				      p,
				      Tag::genTag(my_rank, p, TAG_NODE_GROUP)));
  }

  fillNodeGroupsFromBuffer(communicator, mesh, buffer);

  comm.waitAll(requests);
  comm.freeCommunicationRequest(requests);
  requests.clear();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::synchronizeNodeGroupsSlaves(DistributedSynchronizer & communicator,
							  UInt root,
							  Mesh & mesh) {
  AKANTU_DEBUG_IN();

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  UInt my_rank = comm.whoAmI();


  AKANTU_DEBUG_INFO("Receiving node groups from proc " << root
		    << " TAG("<< Tag::genTag(root, my_rank, TAG_NODE_GROUP) <<")");

  CommunicationStatus status;
  comm.probe<char>(root, Tag::genTag(root, my_rank, TAG_NODE_GROUP), status);

  CommunicationBuffer buffer(status.getSize());
  comm.receive(buffer.storage(), buffer.getSize(), root, Tag::genTag(root, my_rank, TAG_NODE_GROUP));

  fillNodeGroupsFromBuffer(communicator, mesh, buffer);

  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
