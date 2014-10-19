/**
 * @file   filtered_synchronizer.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Mathilde Radiguet <mathilde.radiguet@epfl.ch>
 *
 * @date creation: Wed Sep 18 2013
 * @date last modification: Mon Aug 25 2014
 *
 * @brief  filtered synchronizer
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
FilteredSynchronizer::FilteredSynchronizer(Mesh & mesh,
					   SynchronizerID id,
					   MemoryID memory_id) :
  DistributedSynchronizer(mesh, id, memory_id) {
  AKANTU_DEBUG_IN();
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
FilteredSynchronizer * FilteredSynchronizer::
createFilteredSynchronizer(const DistributedSynchronizer & d_synchronizer,
			   SynchElementFilter & filter) {
  AKANTU_DEBUG_IN();
  
  FilteredSynchronizer & f_synchronizer = 
    *(new FilteredSynchronizer(d_synchronizer.mesh,
			       d_synchronizer.id + ":filtered",
			       d_synchronizer.memory_id));
  
  f_synchronizer.setupSynchronizer(d_synchronizer,
				   filter);

  AKANTU_DEBUG_OUT();
  return &f_synchronizer;
}

/* -------------------------------------------------------------------------- */
void FilteredSynchronizer::setupSynchronizer(const DistributedSynchronizer & d_synchronizer,
					     SynchElementFilter & filter) {
  AKANTU_DEBUG_IN();

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  
  Array<Element> * d_send_element = d_synchronizer.send_element;
  Array<Element> * d_recv_element = d_synchronizer.recv_element;

  std::vector<CommunicationRequest *> isend_requests;

  Array<UInt> ** keep_element_recv = new Array<UInt>*[this->nb_proc];

  // loop over procs
  for (UInt p = 0; p < this->nb_proc; ++p) {
    keep_element_recv[p] = NULL;
    if (p == this->rank) continue;

    // access the element for this proc
    const Array<Element> & unfiltered_elements = d_recv_element[p];
    Array<Element> & filtered_elements = recv_element[p];

    if (unfiltered_elements.getSize() == 0) continue;

    keep_element_recv[p] = new Array<UInt>(0,1,"keep_element_recv");

    // iterator to loop over all source elements
    Array<Element>::const_iterator<Element> it  = unfiltered_elements.begin();
    Array<Element>::const_iterator<Element> end = unfiltered_elements.end();
    
    // if filter accepts this element, push it into the destination elements
    for (UInt el=0; it != end; ++it, ++el) {
      const Element & element = *it;
      if (filter(element)) {
	filtered_elements.push_back(element);
	keep_element_recv[p]->push_back(el);
      }
    }

    keep_element_recv[p]->push_back(-1); // just to be sure to send
				     // something due to some shitty
				     // MPI implementation who do not
				     // know what to do with a 0 size
				     // send

    AKANTU_DEBUG_INFO("I have " << keep_element_recv[p]->getSize() - 1
		      << " elements to still receive from processor " << p
		      << " (communication tag : " 
		      << Tag::genTag(this->rank, 0, RECEIVE_LIST_TAG) << ")");

    isend_requests.push_back(comm.asyncSend<UInt>(keep_element_recv[p]->storage(),
						  keep_element_recv[p]->getSize(),
						  p,
						  Tag::genTag(this->rank, 0, RECEIVE_LIST_TAG)));
  }

  for (UInt p = 0; p < this->nb_proc; ++p) {
    if (p == this->rank) continue;

    const Array<Element> & unfiltered_elements = d_send_element[p];
    Array<Element> & filtered_elements = send_element[p];

    if (unfiltered_elements.getSize() == 0) continue;

    AKANTU_DEBUG_INFO("Waiting list of elements to keep from processor " << p
                      << " (communication tag : " 
		      << Tag::genTag(p, 0, RECEIVE_LIST_TAG) << ")");

    CommunicationStatus status;
    comm.probe<UInt>(p, Tag::genTag(p, 0, RECEIVE_LIST_TAG), status);

    Array<UInt> keep_element(status.getSize(),1,"keep_element_snd");

    AKANTU_DEBUG_INFO("I have " << keep_element.getSize() - 1
                      << " elements to keep in my send list to processor " << p
                      << " (communication tag : " 
		      << Tag::genTag(p, 0, RECEIVE_LIST_TAG) << ")");

    comm.receive(keep_element.storage(),
		 keep_element.getSize(),
		 p,
		 Tag::genTag(p, 0, RECEIVE_LIST_TAG));

    for(UInt i = 0; i < keep_element.getSize() - 1; ++i) {
      filtered_elements.push_back(unfiltered_elements(keep_element(i)));
    }
  }

  comm.waitAll(isend_requests);
  comm.freeCommunicationRequest(isend_requests);

  for (UInt p=0; p < this->nb_proc; ++p) {
    delete keep_element_recv[p];
  }
  delete [] keep_element_recv;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FilteredSynchronizer::updateElementList(Array<Element> * source_elements,
					     Array<Element> * destination_elements,
					     SynchElementFilter & filter) {
  AKANTU_DEBUG_IN();
  
  // loop over procs
  for (UInt p = 0; p < this->nb_proc; ++p) {
    if (p == this->rank) continue;

    // access the element for this proc
    const Array<Element> & unfiltered_elements = source_elements[p];
    Array<Element> & filtered_elements = destination_elements[p];

    // iterator to loop over all source elements
    Array<Element>::const_iterator<Element> it  = unfiltered_elements.begin();
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


__END_AKANTU__
