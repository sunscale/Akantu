/**
 * @file   facet_synchronizer.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Nov 05 2014
 * @date last modification: Fri Jan 26 2018
 *
 * @brief  Facet synchronizer for parallel simulations with cohesive elments
 *
 * @section LICENSE
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
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
        connectivities_for_proc.setID(id + ":connectivities_for_proc:" + std::to_string(proc));
        connectivities_for_proc.initialize(
            mesh, _spatial_dimension = spatial_dimension - 1,
            _with_nb_nodes_per_element = true);

        facet_send_elements.setID(id + ":facet_send_elements:" + std::to_string(proc));
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

} // namespace akantu
