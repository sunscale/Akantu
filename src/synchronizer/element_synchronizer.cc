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
 * @date last modification: Fri Jan 22 2016
 *
 * @brief  implementation of a  communicator using a static_communicator for
 * real
 * send/receive
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
    : SynchronizerImpl<Element>(id, memory_id, mesh.getCommunicator()),
      mesh(mesh), element_to_prank("element_to_prank", id, memory_id) {

  AKANTU_DEBUG_IN();

  if (register_to_event_manager)
    this->mesh.registerEventHandler(*this, event_priority);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
ElementSynchronizer::~ElementSynchronizer() = default;

/* -------------------------------------------------------------------------- */
void ElementSynchronizer::substituteElements(
    const std::map<Element, Element> & old_to_new_elements) {
  auto found_element_end =
    old_to_new_elements.end();

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

  this->removeElements(element_to_remove);
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
  if(not element_to_prank.exists(element.type, element.ghost_type)) {
    // Nicolas: Ok This is nasty I know....
    const_cast<ElementSynchronizer *>(this)->buildElementToPrank();
  }

  return element_to_prank(element);
}

/* -------------------------------------------------------------------------- */
void ElementSynchronizer::reset() {
  AKANTU_DEBUG_IN();

  communications.resetSchemes();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ElementSynchronizer::removeElements(
    const Array<Element> & element_to_remove) {
  AKANTU_DEBUG_IN();

  std::vector<CommunicationRequest> send_requests;
  std::map<UInt, Array<UInt>> list_of_elements_per_proc;

  auto filter_list = [](const Array<UInt> & keep, Array<Element> & list) {
    Array<Element> new_list;
    for (UInt e = 0; e < keep.size() - 1; ++e) {
      Element & el = list(keep(e));
      new_list.push_back(el);
    }
    list.copy(new_list);
  };

  // Handling ghost elements
  for (auto && scheme_pair : communications.iterateSchemes(_recv)) {
    auto proc = scheme_pair.first;
    auto & recv = scheme_pair.second;
    Array<UInt> & keep_list = list_of_elements_per_proc[proc];

    auto rem_it = element_to_remove.begin();
    auto rem_end = element_to_remove.end();

    for (auto && pair : enumerate(recv)) {
      const auto & el = std::get<1>(pair);
      auto pos = std::find(rem_it, rem_end, el);

      if (pos == rem_end) {
        keep_list.push_back(std::get<0>(pair));
      }
    }

    keep_list.push_back(UInt(-1)); // To no send empty arrays

    auto && tag = Tag::genTag(proc, 0, Tag::_MODIFY_SCHEME, this->hash_id);
    AKANTU_DEBUG_INFO("Sending a message of size "
                      << keep_list.size() << " to proc " << proc << " TAG("
                      << tag << ")");
    send_requests.push_back(this->communicator.asyncSend(keep_list, proc, tag));

    auto old_size = recv.size();
    filter_list(keep_list, recv);

    AKANTU_DEBUG_INFO("I had " << old_size << " elements to recv from proc "
                               << proc << " and " << keep_list.size()
                               << " elements to keep. I have " << recv.size()
                               << " elements left.");
  }

  for (auto && scheme_pair : communications.iterateSchemes(_send)) {
    auto proc = scheme_pair.first;
    auto & send = scheme_pair.second;

    CommunicationStatus status;
    auto && tag = Tag::genTag(rank, 0, Tag::_MODIFY_SCHEME, this->hash_id);
    AKANTU_DEBUG_INFO("Getting number of elements of proc "
                      << proc << " to keep - TAG(" << tag << ")");
    this->communicator.probe<UInt>(proc, tag, status);
    Array<UInt> keep_list(status.size());

    AKANTU_DEBUG_INFO("Receiving list of elements ("
                      << keep_list.size() << " elements) to keep for proc "
                      << proc << " TAG(" << tag << ")");
    this->communicator.receive(keep_list, proc, tag);

    auto old_size = send.size();
    filter_list(keep_list, send);

    AKANTU_DEBUG_INFO("I had " << old_size << " elements to send to proc "
                               << proc << " and " << keep_list.size()
                               << " elements to keep. I have " << send.size()
                               << " elements left.");
  }

  this->communicator.waitAll(send_requests);
  this->communicator.freeCommunicationRequest(send_requests);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ElementSynchronizer::renumberElements(
    const ElementTypeMapArray<UInt> & new_numbering) {
  for(auto && sr : iterate_send_recv) {
    for (auto && scheme_pair : communications.iterateSchemes(sr)) {
      auto & list = scheme_pair.second;
      for (auto && el : list) {
        if(new_numbering.exists(el.type, el.ghost_type))
          el.element = new_numbering(el);
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
UInt ElementSynchronizer::sanityCheckDataSize(
    const Array<Element> & elements, const SynchronizationTag &) const {
  return (elements.size() * mesh.getSpatialDimension() * sizeof(Real) +
          sizeof(SynchronizationTag));
}

/* -------------------------------------------------------------------------- */
void ElementSynchronizer::packSanityCheckData(
    CommunicationDescriptor<Element> & comm_desc) const {
  auto & buffer = comm_desc.getBuffer();
  buffer << comm_desc.getTag();

  auto & send_element = comm_desc.getScheme();

  /// pack barycenters in debug mode
  for(auto && element : send_element) {
    Vector<Real> barycenter(mesh.getSpatialDimension());
    mesh.getBarycenter(element.element, element.type, barycenter.storage(),
                       element.ghost_type);
    buffer << barycenter;
  }
}

/* -------------------------------------------------------------------------- */
void ElementSynchronizer::unpackSanityCheckData(
    CommunicationDescriptor<Element> & comm_desc) const {
  auto & buffer = comm_desc.getBuffer();
  const auto & tag = comm_desc.getTag();

  SynchronizationTag t;
  buffer >> t;

  AKANTU_DEBUG_ASSERT(
      t == tag, "The tag received does not correspond to the tag expected");

  auto & recv_element = comm_desc.getScheme();
  auto spatial_dimension = mesh.getSpatialDimension();

  for(auto && element : recv_element) {
    Vector<Real> barycenter_loc(spatial_dimension);
    mesh.getBarycenter(element.element, element.type, barycenter_loc.storage(),
                       element.ghost_type);
    Vector<Real> barycenter(spatial_dimension);
    buffer >> barycenter;
    for (UInt i = 0; i < spatial_dimension; ++i) {
      if (!Math::are_float_equal(barycenter_loc(i), barycenter(i)))
        AKANTU_DEBUG_ERROR("Unpacking an unknown value for the element: "
                           << element << "(barycenter[" << i
                           << "] = " << barycenter_loc(i) << " and buffer[" << i
                           << "] = " << barycenter(i) << ") ["
                           << std::abs(barycenter(i) - barycenter_loc(i))
                           << "] - tag: " << tag);
    }
  }
}

/* -------------------------------------------------------------------------- */
} // namespace akantu
