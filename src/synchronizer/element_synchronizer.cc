/**
 * @file   element_synchronizer.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Dana Christen <dana.christen@epfl.ch>
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Sep 01 2010
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  implementation of a  communicator using a static_communicator for
 * real
 * send/receive
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "element_synchronizer.hh"
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_utils.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
#include <iostream>
#include <map>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
ElementSynchronizer::ElementSynchronizer(Mesh & mesh, const ID & id,
                                         MemoryID memory_id,
                                         bool register_to_event_manager,
                                         EventHandlerPriority event_priority)
    : SynchronizerImpl<Element>(mesh.getCommunicator(), id, memory_id),
      mesh(mesh), element_to_prank("element_to_prank", id, memory_id) {
  AKANTU_DEBUG_IN();

  if (register_to_event_manager)
    this->mesh.registerEventHandler(*this, event_priority);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
ElementSynchronizer::ElementSynchronizer(const ElementSynchronizer & other,
                                         const ID & id,
                                         bool register_to_event_manager,
                                         EventHandlerPriority event_priority)
    : SynchronizerImpl<Element>(other, id), mesh(other.mesh),
      element_to_prank("element_to_prank", id, other.memory_id) {
  AKANTU_DEBUG_IN();

  element_to_prank.copy(other.element_to_prank);

  if (register_to_event_manager)
    this->mesh.registerEventHandler(*this, event_priority);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
ElementSynchronizer::~ElementSynchronizer() = default;

/* -------------------------------------------------------------------------- */
void ElementSynchronizer::substituteElements(
    const std::map<Element, Element> & old_to_new_elements) {
  auto found_element_end = old_to_new_elements.end();

  // substitute old elements with new ones
  for (auto && sr : iterate_send_recv) {
    for (auto && scheme_pair : communications.iterateSchemes(sr)) {
      auto & list = scheme_pair.second;
      for (auto & el : list) {
        auto found_element_it = old_to_new_elements.find(el);
        if (found_element_it != found_element_end)
          el = found_element_it->second;
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void ElementSynchronizer::onElementsChanged(
    const Array<Element> & old_elements_list,
    const Array<Element> & new_elements_list, const ElementTypeMapArray<UInt> &,
    const ChangedElementsEvent &) {
  // create a map to link old elements to new ones
  std::map<Element, Element> old_to_new_elements;

  for (UInt el = 0; el < old_elements_list.size(); ++el) {
    AKANTU_DEBUG_ASSERT(old_to_new_elements.find(old_elements_list(el)) ==
                            old_to_new_elements.end(),
                        "The same element cannot appear twice in the list");

    old_to_new_elements[old_elements_list(el)] = new_elements_list(el);
  }

  substituteElements(old_to_new_elements);

  communications.invalidateSizes();
}

/* -------------------------------------------------------------------------- */
void ElementSynchronizer::onElementsRemoved(
    const Array<Element> & element_to_remove,
    const ElementTypeMapArray<UInt> & new_numbering,
    const RemovedElementsEvent &) {
  AKANTU_DEBUG_IN();

  this->filterScheme([&](auto && element) {
    return std::find(element_to_remove.begin(), element_to_remove.end(),
                     element) == element_to_remove.end();
  });

  this->renumberElements(new_numbering);

  communications.invalidateSizes();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ElementSynchronizer::buildElementToPrank() {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();
  element_to_prank.initialize(mesh, _spatial_dimension = spatial_dimension,
                              _element_kind = _ek_not_defined,
                              _with_nb_element = true, _default_value = rank);

  /// assign prank to all ghost elements
  for (auto && scheme : communications.iterateSchemes(_recv)) {
    auto & recv = scheme.second;
    auto proc = scheme.first;

    for (auto & element : recv) {
      element_to_prank(element) = proc;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Int ElementSynchronizer::getRank(const Element & element) const {
  if (not element_to_prank.exists(element.type, element.ghost_type)) {
    // Nicolas: Ok This is nasty I know....
    const_cast<ElementSynchronizer *>(this)->buildElementToPrank();
  }

  return element_to_prank(element);
}

/* -------------------------------------------------------------------------- */
void ElementSynchronizer::renumberElements(
    const ElementTypeMapArray<UInt> & new_numbering) {
  for (auto && sr : iterate_send_recv) {
    for (auto && scheme_pair : communications.iterateSchemes(sr)) {
      auto & list = scheme_pair.second;
      for (auto && el : list) {
        if (new_numbering.exists(el.type, el.ghost_type))
          el.element = new_numbering(el);
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
UInt ElementSynchronizer::sanityCheckDataSize(const Array<Element> & elements,
                                              const SynchronizationTag & tag,
                                              bool from_comm_desc) const {
  UInt size = SynchronizerImpl<Element>::sanityCheckDataSize(elements, tag,
                                                             from_comm_desc);

  // global connectivities;
  size += mesh.getNbNodesPerElementList(elements) * sizeof(UInt);

  // barycenters
  size += (elements.size() * mesh.getSpatialDimension() * sizeof(Real));
  return size;
}

/* -------------------------------------------------------------------------- */
void ElementSynchronizer::packSanityCheckData(
    CommunicationBuffer & buffer, const Array<Element> & elements,
    const SynchronizationTag & /*tag*/) const {
  for (auto && element : elements) {
    Vector<Real> barycenter(mesh.getSpatialDimension());
    mesh.getBarycenter(element, barycenter);
    buffer << barycenter;

    const auto & conns = mesh.getConnectivity(element.type, element.ghost_type);
    for (auto n : arange(conns.getNbComponent())) {
      buffer << mesh.getNodeGlobalId(conns(element.element, n));
    }
  }
}

/* -------------------------------------------------------------------------- */
void ElementSynchronizer::unpackSanityCheckData(CommunicationBuffer & buffer,
                                                const Array<Element> & elements,
                                                const SynchronizationTag & tag,
                                                UInt proc, UInt rank) const {
  auto spatial_dimension = mesh.getSpatialDimension();

  std::set<SynchronizationTag> skip_conn_tags{
      SynchronizationTag::_smmc_facets_conn,
      SynchronizationTag::_giu_global_conn};

  bool is_skip_tag_conn = skip_conn_tags.find(tag) != skip_conn_tags.end();

  for (auto && element : elements) {
    Vector<Real> barycenter_loc(spatial_dimension);
    mesh.getBarycenter(element, barycenter_loc);

    Vector<Real> barycenter(spatial_dimension);
    buffer >> barycenter;

    auto dist = barycenter_loc.distance(barycenter);
    if (not Math::are_float_equal(dist, 0.)) {
      AKANTU_EXCEPTION("Unpacking an unknown value for the element "
                       << element << "(barycenter " << barycenter_loc
                       << " != buffer " << barycenter << ") [" << dist
                       << "] - tag: " << tag << " comm from " << proc << " to "
                       << rank);
    }

    const auto & conns = mesh.getConnectivity(element.type, element.ghost_type);
    Vector<UInt> global_conn(conns.getNbComponent());
    Vector<UInt> local_global_conn(conns.getNbComponent());

    auto is_same = true;
    for (auto n : arange(global_conn.size())) {
      buffer >> global_conn(n);

      auto node = conns(element.element, n);
      local_global_conn(n) = mesh.getNodeGlobalId(node);

      is_same &= is_skip_tag_conn or mesh.isPureGhostNode(node) or
                 (local_global_conn(n) == global_conn(n));
    }

    if (not is_same) {
      AKANTU_DEBUG_WARNING(
          "The connectivity of the element "
          << element << " " << local_global_conn
          << " does not match the connectivity of the equivalent "
             "element on proc "
          << proc << " " << global_conn << " in communication with tag "
          << tag);
    }
  }
}

/* -------------------------------------------------------------------------- */
} // namespace akantu
