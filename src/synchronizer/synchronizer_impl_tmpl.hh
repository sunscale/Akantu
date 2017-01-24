/**
 * @file   synchronizer_impl_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Aug 24 13:29:47 2016
 *
 * @brief  Implementation of the SynchronizerImpl
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
#include "synchronizer_impl.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SYNCHRONIZER_IMPL_TMPL_HH__
#define __AKANTU_SYNCHRONIZER_IMPL_TMPL_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class Entity>
SynchronizerImpl<Entity>::SynchronizerImpl(const ID & id,
                                           MemoryID memory_id,
                                           const StaticCommunicator & comm)
    : Synchronizer(id, memory_id, comm), communications(comm) {}

/* -------------------------------------------------------------------------- */
template <class Entity>
void SynchronizerImpl<Entity>::asynchronousSynchronizeImpl(
    DataAccessor<Entity> & data_accessor, const SynchronizationTag & tag) {
  AKANTU_DEBUG_IN();

  if (!communications.hasCommunication(tag))
    computeBufferSize(data_accessor, tag);

  // Posting the receive -------------------------------------------------------
  if (communications.hasPendingRecv(tag)) {
    AKANTU_CUSTOM_EXCEPTION_INFO(
        debug::CommunicationException(),
        "There must be some pending receive communications."
            << " Tag is " << tag);
  }

  typedef typename Communications<Entity>::iterator iterator;
  iterator recv_it = communications.begin_recv(tag);
  iterator recv_end = communications.end_recv(tag);
  for (; recv_it != recv_end; ++recv_it) {
    CommunicationDescriptor<Entity> comm_desc = *recv_it;
    comm_desc.postRecv(this->hash_id);
  }

  // Posting the sends -------------------------------------------------------
  if (communications.hasPendingSend(tag)) {
    AKANTU_CUSTOM_EXCEPTION_INFO(
        debug::CommunicationException(),
        "There must be some pending sending communications."
            << " Tag is " << tag);
  }

  iterator send_it = communications.begin_send(tag);
  iterator send_end = communications.end_send(tag);
  for (; send_it != send_end; ++send_it) {
    CommunicationDescriptor<Entity> comm_desc = *send_it;

#ifndef AKANTU_NDEBUG
    this->packSanityCheckData(comm_desc);
#endif

    comm_desc.packData(data_accessor);
    comm_desc.postSend(this->hash_id);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Entity>
void SynchronizerImpl<Entity>::waitEndSynchronizeImpl(
    DataAccessor<Entity> & data_accessor, const SynchronizationTag & tag) {
  AKANTU_DEBUG_IN();

  if (!communications.hasPendingRecv(tag))
    AKANTU_CUSTOM_EXCEPTION_INFO(debug::CommunicationException(),
                                 "No communication with the tag \""
                                     << tag << "\" started");

  typedef typename Communications<Entity>::iterator iterator;
  iterator recv_it;
  iterator recv_end = communications.end_recv(tag);
  while ((recv_it = communications.waitAnyRecv(tag)) != recv_end) {
    CommunicationDescriptor<Entity> comm_desc = *recv_it;
#ifndef AKANTU_NDEBUG
    this->unpackSanityCheckData(comm_desc);
#endif

    comm_desc.unpackData(data_accessor);
    comm_desc.freeRequest();
  }

  communications.waitAllSend(tag);
  communications.freeSendRequests(tag);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Entity>
void SynchronizerImpl<Entity>::computeAllBufferSizes(
    DataAccessor<Entity> & data_accessor) {
  auto it = this->communications.begin_tag();
  auto end = this->communications.end_tag();

  for (; it != end; ++it) {
    SynchronizationTag tag = *it;
    this->computeBufferSize(data_accessor, tag);
  }
}

/* -------------------------------------------------------------------------- */
template <class Entity>
void SynchronizerImpl<Entity>::computeBufferSizeImpl(
    DataAccessor<Entity> & data_accessor, const SynchronizationTag & tag) {
  AKANTU_DEBUG_IN();

  typedef typename Communications<Entity>::scheme_iterator scheme_iterator;
  typedef typename Communications<Entity>::Scheme Scheme;
  scheme_iterator recv_it = communications.begin_recv_scheme();
  scheme_iterator recv_end = communications.end_recv_scheme();
  for (; recv_it != recv_end; ++recv_it) {
    UInt proc = recv_it->first;
    Scheme & scheme = recv_it->second;
    UInt srecv = 0;
#ifndef AKANTU_NDEBUG
    srecv += sanityCheckDataSize(scheme, tag);
#endif
    srecv += data_accessor.getNbData(scheme, tag);
    AKANTU_DEBUG_INFO("I have " << srecv << "(" << printMemorySize<char>(srecv)
                                << " - " << scheme.getSize()
                                << " element(s)) data to receive from " << proc
                                << " for tag " << tag);
    communications.initializeRecvCommunication(tag, proc, srecv);
  }

  scheme_iterator send_it = communications.begin_send_scheme();
  scheme_iterator send_end = communications.end_send_scheme();
  for (; send_it != send_end; ++send_it) {
    UInt proc = send_it->first;
    Scheme & scheme = send_it->second;
    UInt ssend = 0;
#ifndef AKANTU_NDEBUG
    ssend += sanityCheckDataSize(scheme, tag);
#endif
    ssend += data_accessor.getNbData(scheme, tag);
    AKANTU_DEBUG_INFO("I have " << ssend << "(" << printMemorySize<char>(ssend)
                                << " - " << scheme.getSize()
                                << " element(s)) data to send to " << proc
                                << " for tag " << tag);
    communications.initializeSendCommunication(tag, proc, ssend);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Entity>
UInt SynchronizerImpl<Entity>::sanityCheckDataSize(
    __attribute__((unused)) const Array<Entity> & elements,
    __attribute__((unused)) const SynchronizationTag & tag) const {
  return 0;
}

/* -------------------------------------------------------------------------- */
template <class Entity>
void SynchronizerImpl<Entity>::packSanityCheckData(__attribute__((
    unused)) CommunicationDescriptor<Entity> & comm_desc) const {}

/* -------------------------------------------------------------------------- */
template <class Entity>
void SynchronizerImpl<Entity>::unpackSanityCheckData(__attribute__((
    unused)) CommunicationDescriptor<Entity> & comm_desc) const {}

/* -------------------------------------------------------------------------- */
} // akantu

#endif /* __AKANTU_SYNCHRONIZER_IMPL_TMPL_HH__ */
