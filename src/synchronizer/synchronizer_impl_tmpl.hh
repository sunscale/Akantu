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
SynchronizerImpl<Entity>::SynchronizerImpl(const Communicator & comm, const ID & id, MemoryID memory_id)

    : Synchronizer(comm, id, memory_id), communications(comm) {}

/* -------------------------------------------------------------------------- */
template <class Entity>
void SynchronizerImpl<Entity>::synchronizeOnceImpl(
    DataAccessor<Entity> & data_accessor,
    const SynchronizationTag & tag) const {
  // no need to synchronize
  if (this->nb_proc == 1)
    return;

  using CommunicationRequests = std::vector<CommunicationRequest>;
  using CommunicationBuffers = std::map<UInt, CommunicationBuffer>;

  CommunicationRequests send_requests, recv_requests;
  CommunicationBuffers send_buffers, recv_buffers;

  auto postComm = [&](const CommunicationSendRecv & sr,
                      CommunicationBuffers & buffers,
                      CommunicationRequests & requests) -> void {
    for (auto && pair : communications.iterateSchemes(sr)) {
      auto & proc = pair.first;
      const auto & scheme = pair.second;

      auto & buffer = buffers[proc];
      UInt buffer_size = data_accessor.getNbData(scheme, tag);
      buffer.resize(buffer_size);

      if (sr == _recv) {
        requests.push_back(communicator.asyncReceive(
            buffer, proc,
            Tag::genTag(this->rank, 0, Tag::_SYNCHRONIZE, this->hash_id)));
      } else {
        data_accessor.packData(buffer, scheme, tag);
        send_requests.push_back(communicator.asyncSend(
            buffer, proc,
            Tag::genTag(proc, 0, Tag::_SYNCHRONIZE, this->hash_id)));
      }
    }
  };

  // post the receive requests
  postComm(_recv, recv_buffers, recv_requests);

  // post the send data requests
  postComm(_send, send_buffers, send_requests);

  // treat the receive requests
  UInt request_ready;
  while ((request_ready = communicator.waitAny(recv_requests)) != UInt(-1)) {
    CommunicationRequest & req = recv_requests[request_ready];
    UInt proc = req.getSource();

    CommunicationBuffer & buffer = recv_buffers[proc];
    const auto & scheme = this->communications.getScheme(proc, _recv);

    data_accessor.unpackData(buffer, scheme, tag);

    req.free();
    recv_requests.erase(recv_requests.begin() + request_ready);
  }

  communicator.waitAll(send_requests);
  communicator.freeCommunicationRequest(send_requests);
}

/* -------------------------------------------------------------------------- */
template <class Entity>
void SynchronizerImpl<Entity>::asynchronousSynchronizeImpl(
    const DataAccessor<Entity> & data_accessor,
    const SynchronizationTag & tag) {
  AKANTU_DEBUG_IN();

  if (not this->communications.hasCommunicationSize(tag))
    this->computeBufferSize(data_accessor, tag);

  this->communications.incrementCounter(tag);

  // Posting the receive -------------------------------------------------------
  if (this->communications.hasPendingRecv(tag)) {
    AKANTU_CUSTOM_EXCEPTION_INFO(
        debug::CommunicationException(),
        "There must still be some pending receive communications."
            << " Tag is " << tag << " Cannot start new ones");
  }

  for (auto && comm_desc : this->communications.iterateRecv(tag)) {
    comm_desc.postRecv(this->hash_id);
  }

  // Posting the sends -------------------------------------------------------
  if (communications.hasPendingSend(tag)) {
    AKANTU_CUSTOM_EXCEPTION_INFO(
        debug::CommunicationException(),
        "There must be some pending sending communications."
            << " Tag is " << tag);
  }

  for (auto && comm_desc : this->communications.iterateSend(tag)) {
    comm_desc.resetBuffer();

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
  if (this->communications.begin(tag, _recv) !=
      this->communications.end(tag, _recv) &&
      !this->communications.hasPendingRecv(tag))
    AKANTU_CUSTOM_EXCEPTION_INFO(debug::CommunicationException(),
                                 "No pending communication with the tag \""
                                     << tag);
#endif

  auto recv_end = this->communications.end(tag, _recv);
  decltype(recv_end) recv_it;

  while ((recv_it = this->communications.waitAnyRecv(tag)) != recv_end) {
    auto && comm_desc = *recv_it;
#ifndef AKANTU_NDEBUG
    this->unpackSanityCheckData(comm_desc);
#endif

    comm_desc.unpackData(data_accessor);
    comm_desc.resetBuffer();
    comm_desc.freeRequest();
  }

  this->communications.waitAllSend(tag);
  this->communications.freeSendRequests(tag);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Entity>
void SynchronizerImpl<Entity>::computeAllBufferSizes(
    const DataAccessor<Entity> & data_accessor) {
  for (auto && tag : this->communications.iterateTags()) {
    this->computeBufferSize(data_accessor, tag);
  }
}

/* -------------------------------------------------------------------------- */
template <class Entity>
void SynchronizerImpl<Entity>::computeBufferSizeImpl(
    const DataAccessor<Entity> & data_accessor,
    const SynchronizationTag & tag) {
  AKANTU_DEBUG_IN();

  if (not this->communications.hasCommunication(tag)) {
    this->communications.initializeCommunications(tag);
    AKANTU_DEBUG_ASSERT(communications.hasCommunication(tag) == true,
                        "Communications where not properly initialized");
  }

  for (auto sr : iterate_send_recv) {
    for (auto && pair : this->communications.iterateSchemes(sr)) {
      auto proc = pair.first;
      const auto & scheme = pair.second;
      UInt size = 0;
#ifndef AKANTU_NDEBUG
      size += this->sanityCheckDataSize(scheme, tag);
#endif
      size += data_accessor.getNbData(scheme, tag);
      AKANTU_DEBUG_INFO("I have "
                        << size << "(" << printMemorySize<char>(size) << " - "
                        << scheme.size() << " element(s)) data to "
                        << std::string(sr == _recv ? "receive from" : "send to")
                        << proc << " for tag " << tag);

      this->communications.setCommunicationSize(tag, proc, size, sr);
    }
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <typename Entity>
template <typename Pred>
void SynchronizerImpl<Entity>::split(SynchronizerImpl<Entity> & in_synchronizer,
                                     Pred && pred) {
  AKANTU_DEBUG_IN();

  auto filter_list = [&](auto & list, auto & new_list) {
    auto copy = list;
    list.resize(0);
    new_list.resize(0);

    for (auto && entity : copy) {
      if (std::forward<Pred>(pred)(entity)) {
        new_list.push_back(entity);
      } else {
        list.push_back(entity);
      }
    }
  };

  for (auto sr : iterate_send_recv) {
    for (auto & scheme_pair :
         in_synchronizer.communications.iterateSchemes(sr)) {
      auto proc = scheme_pair.first;
      auto & scheme = scheme_pair.second;
      auto & new_scheme = communications.createScheme(proc, sr);
      filter_list(new_scheme, scheme);
    }
  }

  in_synchronizer.communications.invalidateSizes();
  communications.invalidateSizes();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <typename Entity>
template <typename Updater>
void SynchronizerImpl<Entity>::updateSchemes(Updater && scheme_updater) {
  for (auto sr : iterate_send_recv) {
    for (auto & scheme_pair : communications.iterateSchemes(sr)) {
      auto proc = scheme_pair.first;
      auto & scheme = scheme_pair.second;
      std::forward<Updater>(scheme_updater)(scheme, proc, sr);
    }
  }

  communications.invalidateSizes();
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
} // namespace akantu

#endif /* __AKANTU_SYNCHRONIZER_IMPL_TMPL_HH__ */
