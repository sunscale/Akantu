/**
 * @file   filtered_synchronizer.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Mathilde Radiguet <mathilde.radiguet@epfl.ch>
 *
 * @date creation: Wed Sep 18 2013
 * @date last modification: Sat Jul 11 2015
 *
 * @brief  filtered synchronizer
 *
 * @section LICENSE
 *
 * Copyright  (©)  2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de Lausanne)
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
#include "filtered_synchronizer.hh"
/* -------------------------------------------------------------------------- */
#include <map>

/* -------------------------------------------------------------------------- */
namespace akantu {

/* -------------------------------------------------------------------------- */
FilteredSynchronizer::FilteredSynchronizer(Mesh & mesh, SynchronizerID id,
                                           MemoryID memory_id)
    : ElementSynchronizer(mesh, id, memory_id) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
FilteredSynchronizer * FilteredSynchronizer::createFilteredSynchronizer(
    const ElementSynchronizer & d_synchronizer, SynchElementFilter & filter) {
  AKANTU_DEBUG_IN();

  FilteredSynchronizer & f_synchronizer = *(new FilteredSynchronizer(
      d_synchronizer.mesh, d_synchronizer.id + ":filtered",
      d_synchronizer.memory_id));

  f_synchronizer.setupSynchronizer(d_synchronizer, filter);

  AKANTU_DEBUG_OUT();
  return &f_synchronizer;
}

/* -------------------------------------------------------------------------- */
void FilteredSynchronizer::setupSynchronizer(
    const ElementSynchronizer & d_synchronizer, SynchElementFilter & filter) {
  AKANTU_DEBUG_IN();

  std::vector<CommunicationRequest> isend_requests;
  std::map<UInt, Array<UInt> > keep_element_recv;

  auto rs_it = d_synchronizer.communications.begin_recv_scheme();
  auto rs_end = d_synchronizer.communications.end_recv_scheme();

  // Filter the recv list and communicate it to the sender processor
  for (; rs_it != rs_end; ++rs_it) {
    const auto & scheme = rs_it->second;
    auto & proc = rs_it->first;

    Array<UInt> & keep_element = keep_element_recv[proc];
    auto & new_scheme = this->communications.createRecvScheme(proc);

    UInt el_pos = 0;
    for (const auto & element : scheme) {
      if (filter(element)) {
        new_scheme.push_back(element);
        keep_element.push_back(el_pos);
      }
      ++el_pos;
    }

    keep_element.push_back(-1); // just to be sure to send
                                // something due to some shitty
                                // MPI implementation who do not
                                // know what to do with a 0 size
                                // send

    Tag tag = Tag::genTag(this->rank, 0, RECEIVE_LIST_TAG);
    AKANTU_DEBUG_INFO("I have " << keep_element.getSize() - 1
                                << " elements to still receive from processor "
                                << proc << " (communication tag : " << tag
                                << ")");
    isend_requests.push_back(communicator.asyncSend(keep_element, proc, tag));
  }

  auto ss_it = d_synchronizer.communications.begin_send_scheme();
  auto ss_end = d_synchronizer.communications.end_send_scheme();

  // Receive the element to remove from the sender scheme
  for (; ss_it != ss_end; ++ss_it) {
    const auto & scheme = rs_it->second;
    auto & proc = rs_it->first;

    auto & new_scheme = this->communications.createSendScheme(proc);

    Tag tag = Tag::genTag(proc, 0, RECEIVE_LIST_TAG);
    AKANTU_DEBUG_INFO("Waiting list of elements to keep from processor "
                      << proc << " (communication tag : " << tag << ")");

    CommunicationStatus status;
    communicator.probe<UInt>(proc, tag, status);

    Array<UInt> keep_element(status.getSize());

    AKANTU_DEBUG_INFO("I have "
                      << keep_element.getSize() - 1
                      << " elements to keep in my send list to processor "
                      << proc << " (communication tag : " << tag << ")");

    communicator.receive(keep_element, proc, tag);

    for (UInt i = 0; i < keep_element.getSize() - 1; ++i) {
      new_scheme.push_back(scheme(keep_element(i)));
    }
  }

  communicator.waitAll(isend_requests);
  communicator.freeCommunicationRequest(isend_requests);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FilteredSynchronizer::updateElementList(
    Array<Element> * source_elements, Array<Element> * destination_elements,
    SynchElementFilter & filter) {
  AKANTU_DEBUG_IN();

  // loop over procs
  for (UInt p = 0; p < this->nb_proc; ++p) {
    if (p == this->rank)
      continue;

    // access the element for this proc
    const Array<Element> & unfiltered_elements = source_elements[p];
    Array<Element> & filtered_elements = destination_elements[p];

    // iterator to loop over all source elements
    Array<Element>::const_iterator<Element> it = unfiltered_elements.begin();
    Array<Element>::const_iterator<Element> end = unfiltered_elements.end();

    // if filter accepts this element, push it into the destination elements
    for (; it != end; ++it) {
      const Element & element = *it;
      if (filter(element)) {
        filtered_elements.push_back(element);
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

} // akantu
