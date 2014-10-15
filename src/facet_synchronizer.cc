/**
 * @file   facet_synchronizer.cc
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Tue Mar 26 09:55:38 2013
 *
 * @brief  Facet synchronizer for parallel simulations with cohesive elments
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

#include "facet_synchronizer.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
FacetSynchronizer::FacetSynchronizer(Mesh & mesh,
				     SynchronizerID id,
				     MemoryID memory_id) :
  DistributedSynchronizer(mesh, id, memory_id) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
FacetSynchronizer * FacetSynchronizer::
createFacetSynchronizer(DistributedSynchronizer & distributed_synchronizer,
			Mesh & mesh,
			SynchronizerID id,
			MemoryID memory_id) {
  AKANTU_DEBUG_IN();

  FacetSynchronizer & f_synchronizer = *(new FacetSynchronizer(mesh, id, memory_id));

  f_synchronizer.setupFacetSynchronization(distributed_synchronizer);

  AKANTU_DEBUG_OUT();
  return &f_synchronizer;
}

/* -------------------------------------------------------------------------- */
void FacetSynchronizer::updateDistributedSynchronizer(DistributedSynchronizer & distributed_synchronizer,
						      DataAccessor & data_accessor,
						      const Mesh & mesh_cohesive) {
  AKANTU_DEBUG_IN();

  Array<Element> * distrib_send_element = distributed_synchronizer.send_element;
  Array<Element> * distrib_recv_element = distributed_synchronizer.recv_element;

  updateElementList(distrib_send_element, send_element, mesh_cohesive);
  updateElementList(distrib_recv_element, recv_element, mesh_cohesive);

  std::map<SynchronizationTag, Communication>::iterator it
    = distributed_synchronizer.communications.begin();

  std::map<SynchronizationTag, Communication>::iterator end
    = distributed_synchronizer.communications.end();

  for (; it != end; ++it) {
    SynchronizationTag tag = it->first;
    distributed_synchronizer.computeBufferSize(data_accessor, tag);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FacetSynchronizer::updateElementList(Array<Element> * elements,
					  const Array<Element> * facets,
					  const Mesh & mesh_cohesive) {
  AKANTU_DEBUG_IN();

  ElementType current_element_type = _not_defined;
  GhostType current_ghost_type = _casper;

  UInt old_nb_cohesive_elements = 0;
  const Array< std::vector<Element> > * element_to_facet = NULL;

  for (UInt p = 0; p < nb_proc; ++p) {
    const Array<Element> & fa = facets[p];
    Array<Element> & el = elements[p];

    Array<Element>::const_iterator<Element> it = fa.begin();
    Array<Element>::const_iterator<Element> end = fa.end();

    for (; it != end; ++it) {
      const Element & facet = *it;

      if(facet.type != current_element_type || facet.ghost_type != current_ghost_type) {
	current_element_type = facet.type;
	current_ghost_type   = facet.ghost_type;

	ElementType current_coh_element_type
	  = FEEngine::getCohesiveElementType(current_element_type);

	/// compute old number of cohesive elements
	old_nb_cohesive_elements = mesh_cohesive.getNbElement(current_coh_element_type,
							      current_ghost_type);
	old_nb_cohesive_elements -= mesh.getData<UInt>("facet_to_double",
						       current_element_type,
						       current_ghost_type).getSize();

	element_to_facet = & mesh.getData<std::vector<Element> >("element_to_subelement",
								 current_element_type,
								 current_ghost_type);
      }

      const Element & cohesive_element = (*element_to_facet)(facet.element)[1];

      if (cohesive_element.kind == _ek_cohesive &&
	  cohesive_element.element >= old_nb_cohesive_elements)
	el.push_back(cohesive_element);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FacetSynchronizer::setupFacetSynchronization(DistributedSynchronizer & distributed_synchronizer) {
  AKANTU_DEBUG_IN();

  Array<Element> * distrib_send_element = distributed_synchronizer.send_element;
  Array<Element> * distrib_recv_element = distributed_synchronizer.recv_element;

  /// build rank to facet correspondance
  ElementTypeMapArray<UInt> rank_to_facet("rank_to_facet", id);
  initRankToFacet(rank_to_facet);
  buildRankToFacet(rank_to_facet, distrib_recv_element);

  /// generate temp_send/recv element arrays with their connectivity
  Array<ElementTypeMapArray<UInt> *> temp_send_element(nb_proc);
  Array<ElementTypeMapArray<UInt> *> temp_recv_element(nb_proc);
  Array<ElementTypeMapArray<UInt> *> send_connectivity(nb_proc);
  Array<ElementTypeMapArray<UInt> *> recv_connectivity(nb_proc);

  UInt spatial_dimension = mesh.getSpatialDimension();

  for (UInt p = 0; p < nb_proc; ++p) {
    if (p == rank) continue;
    std::stringstream sstr; sstr << p;

    temp_send_element(p) =
      new ElementTypeMapArray<UInt>("temp_send_element_proc_"+sstr.str(), id);
    mesh.initElementTypeMapArray(*temp_send_element(p), 1, spatial_dimension - 1);

    temp_recv_element(p) =
      new ElementTypeMapArray<UInt>("temp_recv_element_proc_"+sstr.str(), id);
    mesh.initElementTypeMapArray(*temp_recv_element(p), 1, spatial_dimension - 1);

    send_connectivity(p) =
      new ElementTypeMapArray<UInt>("send_connectivity_proc_"+sstr.str(), id);
    mesh.initElementTypeMapArray(*send_connectivity(p), 1, spatial_dimension - 1, true);

    recv_connectivity(p) =
      new ElementTypeMapArray<UInt>("recv_connectivity_proc_"+sstr.str(), id);
    mesh.initElementTypeMapArray(*recv_connectivity(p), 1, spatial_dimension - 1, true);
  }

  /// build global connectivity arrays
  getFacetGlobalConnectivity<_not_ghost>(distributed_synchronizer,
					 rank_to_facet,
					 distrib_send_element,
					 send_connectivity,
					 temp_send_element);

  getFacetGlobalConnectivity<_ghost>(distributed_synchronizer,
				     rank_to_facet,
				     distrib_recv_element,
				     recv_connectivity,
				     temp_recv_element);

  /// build send/recv facet arrays
  buildSendElementList(send_connectivity, recv_connectivity, temp_send_element);
  buildRecvElementList(temp_recv_element);

#ifndef AKANTU_NDEBUG
  /// count recv facets for each processor
  Array<UInt> nb_facets_recv(nb_proc);
  nb_facets_recv.clear();

  Mesh::type_iterator first = mesh.firstType(spatial_dimension - 1, _ghost);
  Mesh::type_iterator last  = mesh.lastType(spatial_dimension - 1, _ghost);

  for (; first != last; ++first) {
    const Array<UInt> & r_to_f = rank_to_facet(*first, _ghost);
    UInt nb_facet = r_to_f.getSize();

    for (UInt f = 0; f < nb_facet; ++f) {
      UInt proc = r_to_f(f);
      if (proc != rank)
	++nb_facets_recv(proc);
    }
  }

  for (UInt p = 0; p < nb_proc; ++p) {
    AKANTU_DEBUG_ASSERT(nb_facets_recv(p) == recv_element[p].getSize(),
			"Wrong number of recv facets");
  }

#endif

  for (UInt p = 0; p < nb_proc; ++p) {
    delete temp_send_element(p);
    delete temp_recv_element(p);
    delete send_connectivity(p);
    delete recv_connectivity(p);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FacetSynchronizer::buildSendElementList(const Array<ElementTypeMapArray<UInt> *> & send_connectivity,
					     const Array<ElementTypeMapArray<UInt> *> & recv_connectivity,
					     const Array<ElementTypeMapArray<UInt> *> & temp_send_element) {
  AKANTU_DEBUG_IN();

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();

  UInt spatial_dimension = mesh.getSpatialDimension();

  GhostType ghost_type = _ghost;

  Mesh::type_iterator first = mesh.firstType(spatial_dimension - 1, ghost_type);
  Mesh::type_iterator last  = mesh.lastType(spatial_dimension - 1, ghost_type);

  /// do every communication by element type
  for (; first != last; ++first) {

    ElementType facet_type = *first;

    std::vector<CommunicationRequest *> send_requests;
    UInt * send_size = new UInt[nb_proc];

    /// send asynchronous data
    for (UInt p = 0; p < nb_proc; ++p) {
      if (p == rank) continue;

      const Array<UInt> & recv_conn = (*recv_connectivity(p))(facet_type, _ghost);
      send_size[p] = recv_conn.getSize();

      /// send connectivity size
      send_requests.push_back(comm.asyncSend(send_size + p,
					     1,
					     p,
					     Tag::genTag(rank, p, 0)));

      /// send connectivity data
      send_requests.push_back(comm.asyncSend(recv_conn.storage(),
					     recv_conn.getSize() *
					     recv_conn.getNbComponent(),
					     p,
					     Tag::genTag(rank, p, 1)));
    }

    UInt * recv_size = new UInt[nb_proc];
    UInt nb_nodes_per_facet = Mesh::getNbNodesPerElement(facet_type);

    /// receive data
    for (UInt p = 0; p < nb_proc; ++p) {
      if (p == rank) continue;

      /// receive connectivity size
      comm.receive(recv_size + p, 1, p, Tag::genTag(p, rank, 0));

      Array<UInt> conn_to_match(recv_size[p], nb_nodes_per_facet);

      /// receive connectivity
      comm.receive(conn_to_match.storage(),
		   conn_to_match.getSize() * conn_to_match.getNbComponent(),
		   p,
		   Tag::genTag(p, rank, 1));

      const Array<UInt> & send_conn =
	(*send_connectivity(p))(facet_type, _not_ghost);
      const Array<UInt> & list = (*temp_send_element(p))(facet_type, _not_ghost);
      UInt nb_local_facets = send_conn.getSize();

      AKANTU_DEBUG_ASSERT(nb_local_facets == list.getSize(),
			  "connectivity and facet list have different sizes");

      Array<bool> checked(nb_local_facets);
      checked.clear();

      Element facet(facet_type, 0, _not_ghost, _ek_regular);

      Array<UInt>::iterator<Vector<UInt> > c_to_match_it =
	conn_to_match.begin(nb_nodes_per_facet);
      Array<UInt>::iterator<Vector<UInt> > c_to_match_end =
	conn_to_match.end(nb_nodes_per_facet);

      /// for every sent facet of other processors, find the
      /// corresponding one in the local send connectivity data in
      /// order to build the send_element arrays
      for (; c_to_match_it != c_to_match_end; ++c_to_match_it) {

	Array<UInt>::const_iterator<Vector<UInt> > c_local_it =
	  send_conn.begin(nb_nodes_per_facet);
	Array<UInt>::const_iterator<Vector<UInt> > c_local_end =
	  send_conn.end(nb_nodes_per_facet);

	for (UInt f = 0; f < nb_local_facets; ++f, ++c_local_it) {
	  if (checked(f)) continue;

	  if ( (*c_to_match_it) == (*c_local_it) ) {
	    checked(f) = true;
	    facet.element = list(f);
	    send_element[p].push_back(facet);
	    break;
	  }
	}
	AKANTU_DEBUG_ASSERT(c_local_it != c_local_end, "facet not found");
      }
    }

    /// wait for all communications to be done and free the
    /// communication request array
    comm.waitAll(send_requests);
    comm.freeCommunicationRequest(send_requests);

    delete [] send_size;
    delete [] recv_size;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FacetSynchronizer::buildRecvElementList(const Array<ElementTypeMapArray<UInt> *> & temp_recv_element) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();

  for (UInt p = 0; p < nb_proc; ++p) {
    if (p == rank) continue;

    GhostType ghost_type = _ghost;

    Mesh::type_iterator first = mesh.firstType(spatial_dimension - 1, ghost_type);
    Mesh::type_iterator last  = mesh.lastType(spatial_dimension - 1, ghost_type);

    for (; first != last; ++first) {
      ElementType facet_type = *first;

      const Array<UInt> & list = (*temp_recv_element(p))(facet_type, ghost_type);
      UInt nb_local_facets = list.getSize();

      Element facet(facet_type, 0, ghost_type, _ek_regular);

      for (UInt f = 0; f < nb_local_facets; ++f) {
	facet.element = list(f);
	recv_element[p].push_back(facet);
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FacetSynchronizer::initRankToFacet(ElementTypeMapArray<UInt> & rank_to_facet) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();

  mesh.initElementTypeMapArray(rank_to_facet, 1, spatial_dimension - 1);

  GhostType ghost_type = _ghost;

  Mesh::type_iterator first = mesh.firstType(spatial_dimension - 1, ghost_type);
  Mesh::type_iterator last  = mesh.lastType(spatial_dimension - 1, ghost_type);

  for (; first != last; ++first) {
    ElementType type = *first;
    UInt nb_facet = mesh.getNbElement(type, ghost_type);

    Array<UInt> & rank_to_f = rank_to_facet(type, ghost_type);
    rank_to_f.resize(nb_facet);

    for (UInt f = 0; f < nb_facet; ++f)
      rank_to_f(f) = rank;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FacetSynchronizer::buildRankToFacet(ElementTypeMapArray<UInt> & rank_to_facet,
					 const Array<Element> * elements) {
  AKANTU_DEBUG_IN();

  for (UInt p = 0; p < nb_proc; ++p) {
    if (p == rank) continue;
    const Array<Element> & elem = elements[p];
    UInt nb_element = elem.getSize();

    for (UInt el = 0; el < nb_element; ++el) {
      ElementType type = elem(el).type;
      GhostType gt = elem(el).ghost_type;
      UInt el_index = elem(el).element;

      const Array<Element> & facet_to_element
	= mesh.getSubelementToElement(type, gt);
      UInt nb_facets_per_element = Mesh::getNbFacetsPerElement(type);
      ElementType facet_type = Mesh::getFacetType(type);

      for (UInt f = 0; f < nb_facets_per_element; ++f) {
	const Element & facet = facet_to_element(el_index, f);
	if (facet == ElementNull) continue;
	UInt facet_index = facet.element;
	GhostType facet_gt = facet.ghost_type;

	if (facet_gt == _not_ghost) continue;

	Array<UInt> & t_to_f = rank_to_facet(facet_type, facet_gt);
	if ((p < t_to_f(facet_index)) || (t_to_f(facet_index) == rank))
	  t_to_f(facet_index) = p;
      }
    }
  }

  AKANTU_DEBUG_OUT();
}


__END_AKANTU__
