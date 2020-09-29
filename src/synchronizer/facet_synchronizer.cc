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

#if defined(AKANTU_MODULE)
#define AKANTU_MODULE_SAVE_ AKANTU_MODULE
#undef AKANTU_MODULE
#endif

#define AKANTU_MODULE facet_synchronizer


namespace akantu {

/* -------------------------------------------------------------------------- */
FacetSynchronizer::FacetSynchronizer(
    Mesh & mesh, const ElementSynchronizer & element_synchronizer,
    const ID & id, MemoryID memory_id)
    : ElementSynchronizer(mesh, id, memory_id) {

  auto spatial_dimension = mesh.getSpatialDimension();

  element_to_prank.initialize(mesh, _spatial_dimension = spatial_dimension - 1,
                              _ghost_type = _ghost, _with_nb_element = true,
                              _default_value = rank);

  // Build element to prank
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
        const auto & facet = facets(f);
        if (facet == ElementNull) {
          continue;
        }

        if (facet.ghost_type == _not_ghost) {
          continue;
        }

        auto & facet_rank = element_to_prank(facet);
        if ((proc < UInt(facet_rank)) || (UInt(facet_rank) == rank)) {
          facet_rank = proc;
        }
      }
    }
  }

  ElementTypeMapArray<UInt> facet_global_connectivities(
      "facet_global_connectivities", id, memory_id);
  facet_global_connectivities.initialize(
      mesh, _spatial_dimension = spatial_dimension - 1, _with_nb_element = true,
      _with_nb_nodes_per_element = true);
  mesh.getGlobalConnectivity(facet_global_connectivities);

  // \TODO perhaps a global element numbering might be useful here...
  for (auto type : facet_global_connectivities.elementTypes(_spatial_dimension =
                                                                _all_dimensions,
                   _element_kind = _ek_not_defined, _ghost_type = _not_ghost)) {
    auto & conn = facet_global_connectivities(type, _not_ghost);
    auto conn_view = make_view(conn, conn.getNbComponent());
    std::for_each(conn_view.begin(), conn_view.end(), [&](auto & conn) {
      std::sort(conn.storage(), conn.storage() + conn.size());
    });
  }

  /// init facet check tracking
  ElementTypeMapArray<bool> facet_checked("facet_checked", id, memory_id);
  std::map<UInt, ElementTypeMapArray<UInt>> recv_connectivities;

  /// Generate the recv scheme and connnectivities to send to the other
  /// processors
  for (auto && scheme_pair :
       element_synchronizer.communications.iterateSchemes(_recv)) {
    facet_checked.initialize(mesh, _spatial_dimension = spatial_dimension - 1,
                             _ghost_type = _ghost, _with_nb_element = true,
                             _default_value = false);

    auto proc = scheme_pair.first;
    const auto & elements = scheme_pair.second;
    auto & facet_scheme = communications.createScheme(proc, _recv);

    // this creates empty arrays...
    auto & connectivities_for_proc = recv_connectivities[proc];

    connectivities_for_proc.setID(
        id + ":connectivities_for_proc:" + std::to_string(proc));
    connectivities_for_proc.initialize(
        mesh, _spatial_dimension = spatial_dimension - 1,
        _with_nb_nodes_per_element = true, _ghost_type = _ghost);

    // for every element in the element synchronizer communication scheme,
    // check the facets to see if they should be communicated and create a
    // connectivity array to match with the one other processors might send
    for (auto && element : elements) {
      const auto & facet_to_element =
          mesh.getSubelementToElement(element.type, element.ghost_type);
      Vector<Element> facets = facet_to_element.begin(
          facet_to_element.getNbComponent())[element.element];

      for (UInt f = 0; f < facets.size(); ++f) {
        auto & facet = facets(f);

        // exclude no valid facets
        if (facet == ElementNull) {
          continue;
        }

        // exclude _ghost facet from send scheme and _not_ghost from receive
        if (facet.ghost_type != _ghost) {
          continue;
        }

        // exclude facet from other processors then the one of current
        // interest in case of receive scheme
        if (UInt(element_to_prank(facet)) != proc) {
          continue;
        }

        auto & checked = facet_checked(facet);
        // skip already checked facets
        if (checked) {
          continue;
        }

        checked = true;

        facet_scheme.push_back(facet);

        auto & global_conn =
            facet_global_connectivities(facet.type, facet.ghost_type);
        Vector<UInt> conn =
            global_conn.begin(global_conn.getNbComponent())[facet.element];
        std::sort(conn.storage(), conn.storage() + conn.size());

        connectivities_for_proc(facet.type, facet.ghost_type).push_back(conn);
      }
    }
  }

  std::vector<CommunicationRequest> send_requests;
  /// do every communication by element type
  for (auto && type : mesh.elementTypes(spatial_dimension - 1)) {
    for (auto && pair : recv_connectivities) {
      auto proc = std::get<0>(pair);
      const auto & connectivities_for_proc = std::get<1>(pair);
      auto && tag = Tag::genTag(proc, type, 1337);
      send_requests.push_back(
          communicator.asyncSend(connectivities_for_proc(type, _ghost), proc,
                                 tag, CommunicationMode::_synchronous));
    }

    auto nb_nodes_per_facet = Mesh::getNbNodesPerElement(type);

    communicator.receiveAnyNumber<UInt>(
        send_requests,
        [&](auto && proc, auto && message) {
          auto & local_connectivities =
              facet_global_connectivities(type, _not_ghost);
          auto & send_scheme = communications.createScheme(proc, _send);

          auto conn_view = make_view(local_connectivities, nb_nodes_per_facet);
          auto conn_begin = conn_view.begin();
          auto conn_end = conn_view.end();
          for (const auto & c_to_match :
               make_view(message, nb_nodes_per_facet)) {
            auto it = std::find(conn_begin, conn_end, c_to_match);

            if (it != conn_end) {
              auto facet = Element{type, UInt(it - conn_begin), _not_ghost};
              send_scheme.push_back(facet);
            } else {
              AKANTU_EXCEPTION("No local facet found to send to proc "
                               << proc << " corresponding to " << c_to_match);
            }
          }
        },
        Tag::genTag(rank, type, 1337));
  }
}

} // namespace akantu

#if defined(AKANTU_MODULE_SAVE_)
#undef AKANTU_MODULE
#define AKANTU_MODULE AKANTU_MODULE_SAVE_
#undef AKANTU_MODULE_SAVE_
#endif
