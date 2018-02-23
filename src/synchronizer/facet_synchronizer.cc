/**
 * @file   facet_synchronizer.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 *
 * @brief  Facet synchronizer for parallel simulations with cohesive elments
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */
/* -------------------------------------------------------------------------- */
#include "facet_synchronizer.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
FacetSynchronizer::FacetSynchronizer(
    Mesh & mesh, const ElementSynchronizer & element_synchronizer,
    const ID & id, MemoryID memory_id)
    : ElementSynchronizer(mesh, id, memory_id) {

  const auto & comm = mesh.getCommunicator();
  auto spatial_dimension = mesh.getSpatialDimension();

  element_to_prank.initialize(mesh, _spatial_dimension = spatial_dimension - 1,
                              _ghost_type = _ghost, _with_nb_element = true,
                              _default_value = rank);

  for (auto && scheme_pair :
       element_synchronizer.communications.iterateSchemes(_recv)) {
    auto proc = std::get<0>(scheme_pair);
    const auto & scheme = std::get<1>(scheme_pair);

    for (auto && elem : scheme) {
      const auto & facet_to_element =
          mesh.getSubelementToElement(elem.type, elem.ghost_type);
      Vector<Element> facets = facet_to_element.begin(
          facet_to_element.getNbComponent())[elem.element];

      for (UInt f = 0; f < facets.size(); ++f) {
        auto & facet = facets(f);
        if (facet == ElementNull)
          continue;

        if (facet.ghost_type == _not_ghost)
          continue;

        auto & facet_rank = element_to_prank(facet);
        if ((proc < UInt(facet_rank)) || (UInt(facet_rank) == rank))
          facet_rank = proc;
      }
    }
  }

  ElementTypeMapArray<UInt> facet_global_connectivities(
      "facet_global_connectivities", id, memory_id);
  facet_global_connectivities.initialize(
      mesh, _spatial_dimension = spatial_dimension - 1, _with_nb_element = true,
      _with_nb_nodes_per_element = true);
  mesh.getGlobalConnectivity(facet_global_connectivities);

  /// init facet check tracking
  ElementTypeMapArray<bool> facet_checked("facet_checked", id, memory_id);
  std::vector<CommunicationRequest> send_requests;

  std::map<UInt, ElementTypeMapArray<UInt>> send_facets;
  std::map<UInt, ElementTypeMapArray<UInt>> connectivities;

  for (auto && sr : iterate_send_recv) {
    GhostType interesting_ghost_type = sr == _send ? _not_ghost : _ghost;
    for (auto && scheme_pair :
         element_synchronizer.communications.iterateSchemes(sr)) {
      facet_checked.initialize(mesh, _spatial_dimension = spatial_dimension - 1,
                               _ghost_type = interesting_ghost_type,
                               _with_nb_element = true, _default_value = false);

      auto && proc = scheme_pair.first;
      auto & elements = scheme_pair.second;
      auto & facet_scheme = communications.createScheme(proc, sr);

      // this creates empty arrays...
      auto & facet_send_elements = send_facets[proc];
      auto & connectivities_for_proc = connectivities[proc];
      if (sr == _send) {
        connectivities_for_proc.setID(id + ":connectivities_for_proc:");
        connectivities_for_proc.initialize(
            mesh, _spatial_dimension = spatial_dimension - 1,
            _with_nb_nodes_per_element = true);

        facet_send_elements.setID(id + ":facet_send_elements");
        facet_send_elements.initialize(
            mesh, _spatial_dimension = spatial_dimension - 1);
      }

      for (auto && element : elements) {
        const auto & facet_to_element =
            mesh.getSubelementToElement(element.type, element.ghost_type);
        Vector<Element> facets = facet_to_element.begin(
            facet_to_element.getNbComponent())[element.element];

        for (UInt f = 0; f < facets.size(); ++f) {
          auto & facet = facets(f);
          if (facet == ElementNull)
            continue;

          if (facet.ghost_type != interesting_ghost_type)
            continue;

          if (sr == _recv && UInt(element_to_prank(facet)) != proc)
            continue;

          auto & checked = facet_checked(facet);
          if (checked)
            continue;

          checked = true;

          if (sr == _recv) {
            facet_scheme.push_back(facet);
          } else {
            facet_send_elements(facet.type, facet.ghost_type)
                .push_back(facet.element);
          }

          auto & global_conn =
              facet_global_connectivities(facet.type, facet.ghost_type);
          Vector<UInt> conn =
              global_conn.begin(global_conn.getNbComponent())[facet.element];
          std::sort(conn.storage(), conn.storage() + conn.size());
          connectivities_for_proc(facet.type, facet.ghost_type).push_back(conn);
        }
      }
    }
  }

  /// do every communication by element type
  for (auto && type : mesh.elementTypes(spatial_dimension - 1)) {

    for (auto && pair : connectivities) {
      auto proc = std::get<0>(pair);
      const auto & connectivities_for_proc = std::get<1>(pair);
      auto && tag = Tag::genTag(proc, proc, 1337);
      send_requests.push_back(
          comm.asyncSend(connectivities_for_proc(type, _ghost), proc, tag,
                         CommunicationMode::_synchronous));
    }

    auto nb_nodes_per_facet = Mesh::getNbNodesPerElement(type);

    Array<UInt> buffer;
    Element facet{type, 0, _not_ghost};

    comm.receiveAnyNumber(
        send_requests, buffer,
        [&](auto && proc, auto && message) {
          auto & local_connectivities = connectivities[proc](type, _not_ghost);
          auto & list = send_facets[proc](type, _not_ghost);
          auto & send_scheme = communications.getScheme(proc, _send);

          std::vector<bool> checked(list.size(), false);

          for (auto && c_to_match : make_view(message, nb_nodes_per_facet)) {
#if !defined(AKANTU_NDEBUG)
            bool found = false;
#endif
            for (auto && pair : enumerate(
                     make_view(local_connectivities, nb_nodes_per_facet))) {
              UInt f;
              Vector<UInt> local_conn;
              std::tie(f, local_conn) = pair;

              if (checked[f])
                continue;

              if (c_to_match == local_conn) {
#if !defined(AKANTU_NDEBUG)
                found = true;
#endif
                checked[f] = true;
                facet.element = list(f);
                send_scheme.push_back(facet);
                break;
              }
            }
            AKANTU_DEBUG_ASSERT(found, "facet not found");
          }
        },
        [&]() { return Tag::genTag(rank, rank, 1337); });
  }
}

// /* --------------------------------------------------------------------------
// */ void FacetSynchronizer::setupFacetSynchronization(
//     ElementSynchronizer & distributed_synchronizer) {
//   AKANTU_DEBUG_IN();

//   Array<Element> * distrib_send_element =
//   distributed_synchronizer.send_element; Array<Element> *
//   distrib_recv_element = distributed_synchronizer.recv_element;

//   /// build rank to facet correspondance
//   ElementTypeMapArray<UInt> rank_to_facet("rank_to_facet", id);
//   initRankToFacet(rank_to_facet);
//   buildRankToFacet(rank_to_facet, distrib_recv_element);

//   /// generate temp_send/recv element arrays with their connectivity
//   Array<ElementTypeMapArray<UInt> *> temp_send_element(nb_proc);
//   Array<ElementTypeMapArray<UInt> *> temp_recv_element(nb_proc);
//   Array<ElementTypeMapArray<UInt> *> send_connectivity(nb_proc);
//   Array<ElementTypeMapArray<UInt> *> recv_connectivity(nb_proc);

//   UInt spatial_dimension = mesh.getSpatialDimension();

//   for (UInt p = 0; p < nb_proc; ++p) {
//     if (p == rank)
//       continue;
//     std::stringstream sstr;
//     sstr << p;

//     temp_send_element(p) = new ElementTypeMapArray<UInt>(
//         "temp_send_element_proc_" + sstr.str(), id);
//     mesh.initElementTypeMapArray(*temp_send_element(p), 1,
//                                  spatial_dimension - 1);

//     temp_recv_element(p) = new ElementTypeMapArray<UInt>(
//         "temp_recv_element_proc_" + sstr.str(), id);
//     mesh.initElementTypeMapArray(*temp_recv_element(p), 1,
//                                  spatial_dimension - 1);

//     send_connectivity(p) = new ElementTypeMapArray<UInt>(
//         "send_connectivity_proc_" + sstr.str(), id);
//     mesh.initElementTypeMapArray(*send_connectivity(p), 1,
//                                  spatial_dimension - 1, true);

//     recv_connectivity(p) = new ElementTypeMapArray<UInt>(
//         "recv_connectivity_proc_" + sstr.str(), id);
//     mesh.initElementTypeMapArray(*recv_connectivity(p), 1,
//                                  spatial_dimension - 1, true);
//   }

//   /// build global connectivity arrays
//   getFacetGlobalConnectivity<_not_ghost>(distributed_synchronizer,
//                                          rank_to_facet, distrib_send_element,
//                                          send_connectivity,
//                                          temp_send_element);

//   getFacetGlobalConnectivity<_ghost>(distributed_synchronizer, rank_to_facet,
//                                      distrib_recv_element, recv_connectivity,
//                                      temp_recv_element);

//   /// build send/recv facet arrays
//   buildSendElementList(send_connectivity, recv_connectivity,
//   temp_send_element); buildRecvElementList(temp_recv_element);

// #ifndef AKANTU_NDEBUG
//   /// count recv facets for each processor
//   Array<UInt> nb_facets_recv(nb_proc);
//   nb_facets_recv.clear();

//   Mesh::type_iterator first = mesh.firstType(spatial_dimension - 1, _ghost);
//   Mesh::type_iterator last = mesh.lastType(spatial_dimension - 1, _ghost);

//   for (; first != last; ++first) {
//     const Array<UInt> & r_to_f = rank_to_facet(*first, _ghost);
//     UInt nb_facet = r_to_f.getSize();

//     for (UInt f = 0; f < nb_facet; ++f) {
//       UInt proc = r_to_f(f);
//       if (proc != rank)
//         ++nb_facets_recv(proc);
//     }
//   }

//   for (UInt p = 0; p < nb_proc; ++p) {
//     AKANTU_DEBUG_ASSERT(nb_facets_recv(p) == recv_element[p].getSize(),
//                         "Wrong number of recv facets");
//   }

// #endif

//   for (UInt p = 0; p < nb_proc; ++p) {
//     delete temp_send_element(p);
//     delete temp_recv_element(p);
//     delete send_connectivity(p);
//     delete recv_connectivity(p);
//   }

//   AKANTU_DEBUG_OUT();
// }

// /* --------------------------------------------------------------------------
// */ void FacetSynchronizer::buildSendElementList(
//     const Array<ElementTypeMapArray<UInt> *> & send_connectivity,
//     const Array<ElementTypeMapArray<UInt> *> & recv_connectivity,
//     const Array<ElementTypeMapArray<UInt> *> & temp_send_element) {
//   AKANTU_DEBUG_IN();

//   StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();

//   UInt spatial_dimension = mesh.getSpatialDimension();

//   GhostType ghost_type = _ghost;

//   Mesh::type_iterator first = mesh.firstType(spatial_dimension - 1,
//   ghost_type); Mesh::type_iterator last = mesh.lastType(spatial_dimension -
//   1, ghost_type);

//   /// do every communication by element type
//   for (; first != last; ++first) {

//     ElementType facet_type = *first;

//     std::vector<CommunicationRequest *> send_requests;
//     UInt * send_size = new UInt[nb_proc];

//     /// send asynchronous data
//     for (UInt p = 0; p < nb_proc; ++p) {
//       if (p == rank)
//         continue;

//       const Array<UInt> & recv_conn =
//           (*recv_connectivity(p))(facet_type, _ghost);
//       send_size[p] = recv_conn.getSize();

//       /// send connectivity size
//       send_requests.push_back(
//           comm.asyncSend(send_size + p, 1, p, Tag::genTag(rank, p, 0)));

//       /// send connectivity data
//       send_requests.push_back(comm.asyncSend(
//           recv_conn.storage(), recv_conn.getSize() *
//           recv_conn.getNbComponent(), p, Tag::genTag(rank, p, 1)));
//     }

//     UInt * recv_size = new UInt[nb_proc];
//     UInt nb_nodes_per_facet = Mesh::getNbNodesPerElement(facet_type);

//     /// receive data
//     for (UInt p = 0; p < nb_proc; ++p) {
//       if (p == rank)
//         continue;

//       /// receive connectivity size
//       comm.receive(recv_size + p, 1, p, Tag::genTag(p, rank, 0));

//       Array<UInt> conn_to_match(recv_size[p], nb_nodes_per_facet);

//       /// receive connectivity
//       comm.receive(conn_to_match.storage(),
//                    conn_to_match.getSize() * conn_to_match.getNbComponent(),
//                    p, Tag::genTag(p, rank, 1));

//       const Array<UInt> & send_conn =
//           (*send_connectivity(p))(facet_type, _not_ghost);
//       const Array<UInt> & list =
//           (*temp_send_element(p))(facet_type, _not_ghost);
//       UInt nb_local_facets = send_conn.getSize();

//       AKANTU_DEBUG_ASSERT(nb_local_facets == list.getSize(),
//                           "connectivity and facet list have different
//                           sizes");

//       Array<bool> checked(nb_local_facets);
//       checked.clear();

//       Element facet(facet_type, 0, _not_ghost, _ek_regular);

//       Array<UInt>::iterator<Vector<UInt>> c_to_match_it =
//           conn_to_match.begin(nb_nodes_per_facet);
//       Array<UInt>::iterator<Vector<UInt>> c_to_match_end =
//           conn_to_match.end(nb_nodes_per_facet);

//       /// for every sent facet of other processors, find the
//       /// corresponding one in the local send connectivity data in
//       /// order to build the send_element arrays
//       for (; c_to_match_it != c_to_match_end; ++c_to_match_it) {

//         Array<UInt>::const_iterator<Vector<UInt>> c_local_it =
//             send_conn.begin(nb_nodes_per_facet);
//         Array<UInt>::const_iterator<Vector<UInt>> c_local_end =
//             send_conn.end(nb_nodes_per_facet);

//         for (UInt f = 0; f < nb_local_facets; ++f, ++c_local_it) {
//           if (checked(f))
//             continue;

//           if ((*c_to_match_it) == (*c_local_it)) {
//             checked(f) = true;
//             facet.element = list(f);
//             send_element[p].push_back(facet);
//             break;
//           }
//         }
//         AKANTU_DEBUG_ASSERT(c_local_it != c_local_end, "facet not found");
//       }
//     }

//     /// wait for all communications to be done and free the
//     /// communication request array
//     comm.waitAll(send_requests);
//     comm.freeCommunicationRequest(send_requests);

//     delete[] send_size;
//     delete[] recv_size;
//   }

//   AKANTU_DEBUG_OUT();
// }

// /* --------------------------------------------------------------------------
// */ void FacetSynchronizer::buildRecvElementList(
//     const Array<ElementTypeMapArray<UInt> *> & temp_recv_element) {
//   AKANTU_DEBUG_IN();

//   UInt spatial_dimension = mesh.getSpatialDimension();

//   for (UInt p = 0; p < nb_proc; ++p) {
//     if (p == rank)
//       continue;

//     GhostType ghost_type = _ghost;

//     Mesh::type_iterator first =
//         mesh.firstType(spatial_dimension - 1, ghost_type);
//     Mesh::type_iterator last = mesh.lastType(spatial_dimension - 1,
//     ghost_type);

//     for (; first != last; ++first) {
//       ElementType facet_type = *first;

//       const Array<UInt> & list =
//           (*temp_recv_element(p))(facet_type, ghost_type);
//       UInt nb_local_facets = list.getSize();

//       Element facet(facet_type, 0, ghost_type, _ek_regular);

//       for (UInt f = 0; f < nb_local_facets; ++f) {
//         facet.element = list(f);
//         recv_element[p].push_back(facet);
//       }
//     }
//   }

//   AKANTU_DEBUG_OUT();
// }

// void FacetSynchronizer::initRankToFacet(
//     ElementTypeMapArray<UInt> & rank_to_facet) {
//   AKANTU_DEBUG_IN();

//   auto spatial_dimension = mesh.getSpatialDimension();

//   for (; first != last; ++first) {
//     ElementType type = *first;
//     UInt nb_facet = mesh.getNbElement(type, ghost_type);

//     Array<UInt> & rank_to_f = rank_to_facet(type, ghost_type);
//     rank_to_f.resize(nb_facet);

//     for (UInt f = 0; f < nb_facet; ++f)
//       rank_to_f(f) = rank;
//   }

//   AKANTU_DEBUG_OUT();
// }

// /* --------------------------------------------------------------------------
// */ void FacetSynchronizer::buildRankToFacet(
//     ElementTypeMapArray<UInt> & rank_to_facet,
//     const Array<Element> * elements) {
//   AKANTU_DEBUG_IN();

//   for (UInt p = 0; p < nb_proc; ++p) {
//     if (p == rank)
//       continue;
//     const Array<Element> & elem = elements[p];
//     UInt nb_element = elem.getSize();

//     for (UInt el = 0; el < nb_element; ++el) {
//       ElementType type = elem(el).type;
//       GhostType gt = elem(el).ghost_type;
//       UInt el_index = elem(el).element;

//       const Array<Element> & facet_to_element =
//           mesh.getSubelementToElement(type, gt);
//       UInt nb_facets_per_element = Mesh::getNbFacetsPerElement(type);
//       ElementType facet_type = Mesh::getFacetType(type);

//       for (UInt f = 0; f < nb_facets_per_element; ++f) {
//         const Element & facet = facet_to_element(el_index, f);
//         if (facet == ElementNull)
//           continue;
//         UInt facet_index = facet.element;
//         GhostType facet_gt = facet.ghost_type;

//         if (facet_gt == _not_ghost)
//           continue;

//         Array<UInt> & t_to_f = rank_to_facet(facet_type, facet_gt);
//         if ((p < t_to_f(facet_index)) || (t_to_f(facet_index) == rank))
//           t_to_f(facet_index) = p;
//       }
//     }
//   }

//   AKANTU_DEBUG_OUT();
// }

} // namespace akantu
