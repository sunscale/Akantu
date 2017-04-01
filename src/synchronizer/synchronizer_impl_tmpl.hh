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
SynchronizerImpl<Entity>::SynchronizerImpl(const ID & id, MemoryID memory_id,
                                           const StaticCommunicator & comm)
    : Synchronizer(id, memory_id, comm), communications(comm) {}

/* -------------------------------------------------------------------------- */
template <class Entity>
void SynchronizerImpl<Entity>::asynchronousSynchronizeImpl(
    DataAccessor<Entity> & data_accessor, const SynchronizationTag & tag) {
  AKANTU_DEBUG_IN();

  if (!this->communications.hasCommunication(tag))
    this->computeBufferSize(data_accessor, tag);

  this->communications.incrementCounter(tag);

  // Posting the receive -------------------------------------------------------
  if (this->communications.hasPendingRecv(tag)) {
    AKANTU_CUSTOM_EXCEPTION_INFO(
        debug::CommunicationException(),
        "There must be some pending receive communications."
            << " Tag is " << tag);
  }

  auto recv_it = this->communications.begin_recv(tag);
  auto recv_end = this->communications.end_recv(tag);
  for (; recv_it != recv_end; ++recv_it) {
    auto comm_desc = *recv_it;
    comm_desc.postRecv(this->hash_id);
  }

  // Posting the sends -------------------------------------------------------
  if (communications.hasPendingSend(tag)) {
    AKANTU_CUSTOM_EXCEPTION_INFO(
        debug::CommunicationException(),
        "There must be some pending sending communications."
            << " Tag is " << tag);
  }

  auto send_it = communications.begin_send(tag);
  auto send_end = communications.end_send(tag);
  for (; send_it != send_end; ++send_it) {
    auto comm_desc = *send_it;

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

#ifndef AKANTU_NDEBUG
  if (!this->communications.hasPendingRecv(tag) && nb_proc > 1)
    AKANTU_CUSTOM_EXCEPTION_INFO(debug::CommunicationException(),
                                 "No pending communication with the tag \""
                                     << tag);
#endif

  auto recv_end = this->communications.end_recv(tag);
  decltype(recv_end) recv_it;

  while ((recv_it = this->communications.waitAnyRecv(tag)) != recv_end) {
    auto comm_desc = *recv_it;
#ifndef AKANTU_NDEBUG
    this->unpackSanityCheckData(comm_desc);
#endif

    comm_desc.unpackData(data_accessor);
    comm_desc.freeRequest();
  }

  this->communications.waitAllSend(tag);
  this->communications.freeSendRequests(tag);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Entity>
void SynchronizerImpl<Entity>::computeAllBufferSizes(
    DataAccessor<Entity> & data_accessor) {
  auto it = this->communications.begin_tag();
  auto end = this->communications.end_tag();

  for (; it != end; ++it) {
    auto tag = *it;
    this->computeBufferSize(data_accessor, tag);
  }
}

/* -------------------------------------------------------------------------- */
template <class Entity>
void SynchronizerImpl<Entity>::computeBufferSizeImpl(
    DataAccessor<Entity> & data_accessor, const SynchronizationTag & tag) {
  AKANTU_DEBUG_IN();

  this->communications.initializeCommunications(tag);
  AKANTU_DEBUG_ASSERT(communications.hasCommunication(tag) == true,
                      "Communications where not properly initialized");

  auto recv_it = this->communications.begin_recv_scheme();
  auto recv_end = this->communications.end_recv_scheme();
  for (; recv_it != recv_end; ++recv_it) {
    auto proc = recv_it->first;
    auto & scheme = recv_it->second;
    UInt srecv = 0;
#ifndef AKANTU_NDEBUG
    srecv += this->sanityCheckDataSize(scheme, tag);
#endif
    srecv += data_accessor.getNbData(scheme, tag);
    AKANTU_DEBUG_INFO("I have " << srecv << "(" << printMemorySize<char>(srecv)
                                << " - " << scheme.getSize()
                                << " element(s)) data to receive from " << proc
                                << " for tag " << tag);
    this->communications.setRecvCommunicationSize(tag, proc, srecv);
  }

  auto send_it = this->communications.begin_send_scheme();
  auto send_end = this->communications.end_send_scheme();
  for (; send_it != send_end; ++send_it) {
    auto proc = send_it->first;
    auto & scheme = send_it->second;
    UInt ssend = 0;
#ifndef AKANTU_NDEBUG
    ssend += this->sanityCheckDataSize(scheme, tag);
#endif
    ssend += data_accessor.getNbData(scheme, tag);
    AKANTU_DEBUG_INFO("I have " << ssend << "(" << printMemorySize<char>(ssend)
                                << " - " << scheme.getSize()
                                << " element(s)) data to send to " << proc
                                << " for tag " << tag);
    this->communications.setSendCommunicationSize(tag, proc, ssend);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Entity>
UInt SynchronizerImpl<Entity>::sanityCheckDataSize(
    const Array<Entity> &, const SynchronizationTag &) const {
  return 0;
}

/* -------------------------------------------------------------------------- */
template <class Entity>
void SynchronizerImpl<Entity>::packSanityCheckData(
    CommunicationDescriptor<Entity> &) const {}

/* -------------------------------------------------------------------------- */
template <class Entity>
void SynchronizerImpl<Entity>::unpackSanityCheckData(
    CommunicationDescriptor<Entity> &) const {}

/* -------------------------------------------------------------------------- */
} // akantu

#endif /* __AKANTU_SYNCHRONIZER_IMPL_TMPL_HH__ */
