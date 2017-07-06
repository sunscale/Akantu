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

#if defined(AKANTU_DEBUG_TOOLS)
#include "aka_debug_tools.hh"
#endif

/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
ElementSynchronizer::ElementSynchronizer(Mesh & mesh,
                                         const StaticCommunicator & comm,
                                         const ID & id, MemoryID memory_id,
                                         bool register_to_event_manager,
                                         UInt event_priority)
    : SynchronizerImpl<Element>(id, memory_id, comm), mesh(mesh),
      prank_to_element("prank_to_element", id, memory_id) {

  AKANTU_DEBUG_IN();

  if (register_to_event_manager)
    this->mesh.registerEventHandler(*this, event_priority);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
ElementSynchronizer::ElementSynchronizer(Mesh & mesh,
                                         const ID & id, MemoryID memory_id,
                                         bool register_to_event_manager,
                                         UInt event_priority)
    : ElementSynchronizer(mesh, mesh.getCommunicator(), id, memory_id,
                          register_to_event_manager, event_priority) {}

/* -------------------------------------------------------------------------- */
ElementSynchronizer::~ElementSynchronizer() = default;

/* -------------------------------------------------------------------------- */
void ElementSynchronizer::substituteElements(
    const std::map<Element, Element> & old_to_new_elements) {
  // substitute old elements with new ones

  auto subsitute =
      [&old_to_new_elements](Communications<Element>::scheme_iterator it,
                             Communications<Element>::scheme_iterator end) {
        std::map<Element, Element>::const_iterator found_element_it;
        std::map<Element, Element>::const_iterator found_element_end =
            old_to_new_elements.end();

        for (; it != end; ++it) {
          Array<Element> & list = it->second;
          for (UInt el = 0; el < list.getSize(); ++el) {
            found_element_it = old_to_new_elements.find(list(el));
            if (found_element_it != found_element_end)
              list(el) = found_element_it->second;
          }
        }
      };

  subsitute(communications.begin_recv_scheme(),
            communications.end_recv_scheme());
  subsitute(communications.begin_send_scheme(),
            communications.end_send_scheme());
}

/* -------------------------------------------------------------------------- */
void ElementSynchronizer::onElementsChanged(
    const Array<Element> & old_elements_list,
    const Array<Element> & new_elements_list,
    __attribute__((unused)) const ElementTypeMapArray<UInt> & new_numbering,
    __attribute__((unused)) const ChangedElementsEvent & event) {
  // create a map to link old elements to new ones
  std::map<Element, Element> old_to_new_elements;

  for (UInt el = 0; el < old_elements_list.getSize(); ++el) {
    AKANTU_DEBUG_ASSERT(old_to_new_elements.find(old_elements_list(el)) ==
                            old_to_new_elements.end(),
                        "The same element cannot appear twice in the list");

    old_to_new_elements[old_elements_list(el)] = new_elements_list(el);
  }

  substituteElements(old_to_new_elements);
}

/* -------------------------------------------------------------------------- */
void ElementSynchronizer::onElementsRemoved(
    const Array<Element> & element_to_remove,
    const ElementTypeMapArray<UInt> & new_numbering,
    __attribute__((unused)) const RemovedElementsEvent & event) {
  AKANTU_DEBUG_IN();
  this->removeElements(element_to_remove);
  this->renumberElements(new_numbering);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ElementSynchronizer::buildPrankToElement() {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();
  prank_to_element.initialize(mesh, _spatial_dimension = spatial_dimension,
                              _element_kind = _ek_not_defined,
                              _with_nb_element = true,
                              _default_value = rank);

  /// assign prank to all ghost elements
  Communications<Element>::scheme_iterator recv_it =
      communications.begin_recv_scheme();
  Communications<Element>::scheme_iterator recv_end =
      communications.end_recv_scheme();
  for (; recv_it != recv_end; ++recv_it) {
    auto & recv = recv_it->second;
    auto proc = recv_it->first;

    for (auto & element : recv) {
      auto & prank_to_el = prank_to_element(element.type, element.ghost_type);
      prank_to_el(element.element) = proc;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ElementSynchronizer::filterElementsByKind(
    ElementSynchronizer * new_synchronizer, ElementKind kind) {
  AKANTU_DEBUG_IN();

  auto filter_list = [&kind](Array<Element> & list, Array<Element> & new_list) {
    list.resize(0);
    new_list.resize(0);

    Array<Element> copy = list;
    Array<Element>::const_scalar_iterator it = copy.begin();
    Array<Element>::const_scalar_iterator end = copy.end();
    for (; it != end; ++it) {
      const Element & el = *it;
      if (el.kind == kind) {
        new_list.push_back(el);
      } else {
        list.push_back(el);
      }
    }
  };

  Communications<Element>::scheme_iterator recv_it =
      communications.begin_recv_scheme();
  Communications<Element>::scheme_iterator recv_end =
      communications.end_recv_scheme();
  for (; recv_it != recv_end; ++recv_it) {
    UInt proc = recv_it->first;
    Array<Element> & recv = recv_it->second;
    Array<Element> & new_recv =
        new_synchronizer->communications.createRecvScheme(proc);

    filter_list(recv, new_recv);
  }

  Communications<Element>::scheme_iterator send_it =
      communications.begin_send_scheme();
  Communications<Element>::scheme_iterator send_end =
      communications.end_send_scheme();
  for (; send_it != send_end; ++send_it) {
    UInt proc = send_it->first;
    Array<Element> & send = send_it->second;
    Array<Element> & new_send =
        new_synchronizer->communications.createSendScheme(proc);

    filter_list(send, new_send);
  }

  AKANTU_DEBUG_OUT();
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
    for (UInt e = 0; e < keep.getSize() - 1; ++e) {
      Element & el = list(keep(e));
      new_list.push_back(el);
    }
    list.copy(new_list);
  };

  // Handling ghost elements
  Communications<Element>::scheme_iterator recv_it =
      communications.begin_recv_scheme();
  Communications<Element>::scheme_iterator recv_end =
      communications.end_recv_scheme();
  for (; recv_it != recv_end; ++recv_it) {
    Array<Element> & recv = recv_it->second;
    UInt proc = recv_it->first;

    Array<UInt> & keep_list = list_of_elements_per_proc[proc];

    Array<Element>::const_iterator<Element> rem_it = element_to_remove.begin();
    Array<Element>::const_iterator<Element> rem_end = element_to_remove.end();

    Array<Element>::const_scalar_iterator it = recv.begin();
    Array<Element>::const_scalar_iterator end = recv.end();
    for (UInt e = 0; it != end; ++it, ++e) {
      const Element & el = *it;
      Array<Element>::const_iterator<Element> pos =
          std::find(rem_it, rem_end, el);
      if (pos == rem_end) {
        keep_list.push_back(e);
      }
    }

    keep_list.push_back(UInt(-1)); // To no send empty arrays

    Tag tag = Tag::genTag(proc, 0, Tag::_MODIFY_SCHEME, this->hash_id);
    AKANTU_DEBUG_INFO("Sending a message of size "
                      << keep_list.getSize() << " to proc " << proc << " TAG("
                      << tag << ")");
    send_requests.push_back(this->communicator.asyncSend(keep_list, proc, tag));

    UInt old_size = recv.getSize();
    filter_list(keep_list, recv);

    AKANTU_DEBUG_INFO("I had " << old_size << " elements to recv from proc "
                               << proc << " and " << keep_list.getSize()
                               << " elements to keep. I have " << recv.getSize()
                               << " elements left.");
  }

  Communications<Element>::scheme_iterator send_it =
      communications.begin_send_scheme();
  Communications<Element>::scheme_iterator send_end =
      communications.end_send_scheme();
  for (; send_it != send_end; ++send_it) {
    UInt proc = send_it->first;
    Array<Element> send = send_it->second;

    CommunicationStatus status;
    Tag tag = Tag::genTag(rank, 0, Tag::_MODIFY_SCHEME, this->hash_id);
    AKANTU_DEBUG_INFO("Getting number of elements of proc "
                      << proc << " to keep - TAG(" << tag << ")");
    this->communicator.probe<UInt>(proc, tag, status);
    Array<UInt> keep_list(status.getSize());

    AKANTU_DEBUG_INFO("Receiving list of elements ("
                      << keep_list.getSize() << " elements) to keep for proc "
                      << proc << " TAG(" << tag << ")");
    this->communicator.receive(keep_list, proc, tag);

    UInt old_size = send.getSize();
    filter_list(keep_list, send);

    AKANTU_DEBUG_INFO("I had " << old_size << " elements to send to proc "
                               << proc << " and " << keep_list.getSize()
                               << " elements to keep. I have " << send.getSize()
                               << " elements left.");
  }

  this->communicator.waitAll(send_requests);
  this->communicator.freeCommunicationRequest(send_requests);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ElementSynchronizer::renumberElements(
    const ElementTypeMapArray<UInt> & new_numbering) {

  auto renumber =
      [&new_numbering](Communications<Element>::scheme_iterator it,
                       Communications<Element>::scheme_iterator end) {
        for (; it != end; ++it) {
          Array<Element> & list = it->second;
          Array<Element>::scalar_iterator el_it = list.begin();
          Array<Element>::scalar_iterator el_end = list.end();
          for (; el_it != el_end; ++el_it) {
            Element & el = *el_it;
            el.element = new_numbering(el.type, el.ghost_type)(el.element);
          }
        }
      };

  renumber(communications.begin_recv_scheme(),
           communications.end_recv_scheme());
  renumber(communications.begin_send_scheme(),
           communications.end_send_scheme());
}

/* -------------------------------------------------------------------------- */
UInt ElementSynchronizer::sanityCheckDataSize(
    const Array<Element> & elements, const SynchronizationTag &) const {
  return (elements.getSize() * mesh.getSpatialDimension() * sizeof(Real) +
          sizeof(SynchronizationTag));
}

/* -------------------------------------------------------------------------- */
void ElementSynchronizer::packSanityCheckData(
    CommunicationDescriptor<Element> & comm_desc) const {
  CommunicationBuffer & buffer = comm_desc.getBuffer();
  buffer << comm_desc.getTag();

  Communications<Element>::Scheme & send_element = comm_desc.getScheme();

  /// pack barycenters in debug mode
  Array<Element>::const_iterator<Element> bit = send_element.begin();
  Array<Element>::const_iterator<Element> bend = send_element.end();
  for (; bit != bend; ++bit) {
    const Element & element = *bit;
    Vector<Real> barycenter(mesh.getSpatialDimension());
    mesh.getBarycenter(element.element, element.type, barycenter.storage(),
                       element.ghost_type);
    buffer << barycenter;
  }
}

/* -------------------------------------------------------------------------- */
void ElementSynchronizer::unpackSanityCheckData(
    CommunicationDescriptor<Element> & comm_desc) const {
  CommunicationBuffer & buffer = comm_desc.getBuffer();
  const SynchronizationTag & tag = comm_desc.getTag();

  SynchronizationTag t;
  buffer >> t;

  AKANTU_DEBUG_ASSERT(
      t == tag, "The tag received does not correspond to the tag expected");

  Communications<Element>::Scheme & recv_element = comm_desc.getScheme();
  Array<Element>::const_iterator<Element> bit = recv_element.begin();
  Array<Element>::const_iterator<Element> bend = recv_element.end();

  UInt spatial_dimension = mesh.getSpatialDimension();

  for (; bit != bend; ++bit) {
    const Element & element = *bit;

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
