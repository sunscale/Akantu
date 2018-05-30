/**
 * @file   synchronizer_impl_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Sep 07 2016
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Implementation of the SynchronizerImpl
 *
 * @section LICENSE
 *
 * Copyright (©) 2016-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "synchronizer_impl.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SYNCHRONIZER_IMPL_TMPL_HH__
#define __AKANTU_SYNCHRONIZER_IMPL_TMPL_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class Entity>
SynchronizerImpl<Entity>::SynchronizerImpl(const Communicator & comm,
                                           const ID & id, MemoryID memory_id)

    : Synchronizer(comm, id, memory_id), communications(comm) {}

/* -------------------------------------------------------------------------- */
template <class Entity>
SynchronizerImpl<Entity>::SynchronizerImpl(const SynchronizerImpl & other,
                                           const ID & id)
    : Synchronizer(other), communications(other.communications) {
  this->id = id;
}

/* -------------------------------------------------------------------------- */
template <class Entity>
void SynchronizerImpl<Entity>::communicateOnce(
    const std::tuple<CommunicationSendRecv, CommunicationSendRecv> &
        send_recv_schemes,
    const Tag::CommTags & comm_tag, DataAccessor<Entity> & data_accessor,
    const SynchronizationTag & tag) const {
  // no need to synchronize
  if (this->nb_proc == 1)
    return;

  CommunicationSendRecv send_dir, recv_dir;
  std::tie(send_dir, recv_dir) = send_recv_schemes;

  using CommunicationRequests = std::vector<CommunicationRequest>;
  using CommunicationBuffers = std::map<UInt, CommunicationBuffer>;

  CommunicationRequests send_requests, recv_requests;
  CommunicationBuffers send_buffers, recv_buffers;

  auto postComm = [&](const auto & sr, auto & buffers,
                      auto & requests) -> void {
    for (auto && pair : communications.iterateSchemes(sr)) {
      auto & proc = pair.first;
      const auto & scheme = pair.second;

      if (scheme.size() == 0)
        continue;

      auto & buffer = buffers[proc];

      auto buffer_size = data_accessor.getNbData(scheme, tag);
      if (buffer_size == 0)
        continue;

#ifndef AKANTU_NDEBUG
      buffer_size += this->sanityCheckDataSize(scheme, tag, false);
#endif

      buffer.resize(buffer_size);

      if (sr == recv_dir) {
        requests.push_back(communicator.asyncReceive(
            buffer, proc,
            Tag::genTag(this->rank, tag, comm_tag, this->hash_id)));
      } else {
#ifndef AKANTU_NDEBUG
        this->packSanityCheckData(buffer, scheme, tag);
#endif
        data_accessor.packData(buffer, scheme, tag);

        AKANTU_DEBUG_ASSERT(
            buffer.getPackedSize() == buffer.size(),
            "The data accessor did not pack all the data it "
            "promised  in communication with tag "
                << tag << " (Promised: " << buffer.size()
                << "bytes, packed: " << buffer.getPackedSize() << "bytes [avg: "
                << Real(buffer.size() - buffer.getPackedSize()) / scheme.size()
                << "bytes per entity missing])");

        send_requests.push_back(communicator.asyncSend(
            buffer, proc, Tag::genTag(proc, tag, comm_tag, this->hash_id)));
      }
    }
  };

  // post the receive requests
  postComm(recv_dir, recv_buffers, recv_requests);

  // post the send data requests
  postComm(send_dir, send_buffers, send_requests);

  // treat the receive requests
  UInt request_ready;
  while ((request_ready = communicator.waitAny(recv_requests)) != UInt(-1)) {
    auto & req = recv_requests[request_ready];
    auto proc = req.getSource();

    auto & buffer = recv_buffers[proc];
    const auto & scheme = this->communications.getScheme(proc, recv_dir);

#ifndef AKANTU_NDEBUG
    this->unpackSanityCheckData(buffer, scheme, tag, proc, this->rank);
#endif

    data_accessor.unpackData(buffer, scheme, tag);

    AKANTU_DEBUG_ASSERT(
        buffer.getLeftToUnpack() == 0,
        "The data accessor ignored some data in communication with tag "
            << tag);

    req.free();
    recv_requests.erase(recv_requests.begin() + request_ready);
  }

  communicator.waitAll(send_requests);
  communicator.freeCommunicationRequest(send_requests);
}

/* -------------------------------------------------------------------------- */
template <class Entity>
void SynchronizerImpl<Entity>::slaveReductionOnceImpl(
    DataAccessor<Entity> & data_accessor,
    const SynchronizationTag & tag) const {
  communicateOnce(std::make_tuple(_recv, _send), Tag::_REDUCE, data_accessor,
                  tag);
}

/* -------------------------------------------------------------------------- */
template <class Entity>
void SynchronizerImpl<Entity>::synchronizeOnceImpl(
    DataAccessor<Entity> & data_accessor,
    const SynchronizationTag & tag) const {
  communicateOnce(std::make_tuple(_send, _recv), Tag::_SYNCHRONIZE,
                  data_accessor, tag);
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
template <typename Entity> void SynchronizerImpl<Entity>::reset() {
  AKANTU_DEBUG_IN();
  communications.resetSchemes();
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
      filter_list(scheme, new_scheme);
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
template <typename Entity>
template <typename Pred>
void SynchronizerImpl<Entity>::filterScheme(Pred && pred) {
  std::vector<CommunicationRequest> requests;
  std::unordered_map<UInt, Array<UInt>> keep_entities;

  auto filter_list = [](const auto & keep, auto & list) {
    Array<Entity> new_list;
    for (const auto & keep_entity : keep) {
      const Entity & entity = list(keep_entity);
      new_list.push_back(entity);
    }
    list.copy(new_list);
  };

  // loop over send_schemes
  for (auto & scheme_pair : communications.iterateSchemes(_recv)) {
    auto proc = scheme_pair.first;
    auto & scheme = scheme_pair.second;

    auto & keep_entity = keep_entities[proc];
    for (auto && entity : enumerate(scheme)) {
      if (pred(std::get<1>(entity))) {
        keep_entity.push_back(std::get<0>(entity));
      }
    }

    auto tag = Tag::genTag(this->rank, 0, Tag::_MODIFY_SCHEME);
    AKANTU_DEBUG_INFO("I have " << keep_entity.size()
                                << " elements to still receive from processor "
                                << proc << " (communication tag : " << tag
                                << ")");

    filter_list(keep_entity, scheme);
    requests.push_back(communicator.asyncSend(keep_entity, proc, tag));
  }

  // clean the receive scheme
  for (auto & scheme_pair : communications.iterateSchemes(_send)) {
    auto proc = scheme_pair.first;
    auto & scheme = scheme_pair.second;

    auto tag = Tag::genTag(proc, 0, Tag::_MODIFY_SCHEME);
    AKANTU_DEBUG_INFO("Waiting list of elements to keep from processor "
                      << proc << " (communication tag : " << tag << ")");

    CommunicationStatus status;
    communicator.probe<UInt>(proc, tag, status);

    Array<UInt> keep_entity(status.size(), 1, "keep_element");
    AKANTU_DEBUG_INFO("I have "
                      << keep_entity.size()
                      << " elements to keep in my send list to processor "
                      << proc << " (communication tag : " << tag << ")");

    communicator.receive(keep_entity, proc, tag);

    filter_list(keep_entity, scheme);
  }

  communicator.waitAll(requests);
  communicator.freeCommunicationRequest(requests);
  communications.invalidateSizes();
}

/* -------------------------------------------------------------------------- */
template <class Entity> void SynchronizerImpl<Entity>::swapSendRecv() {
  communications.swapSendRecv();
}

/* -------------------------------------------------------------------------- */
template <class Entity>
void SynchronizerImpl<Entity>::
copySchemes(const SynchronizerImpl & other) {
  reset();

  for (auto sr : iterate_send_recv) {
    for (auto & scheme_pair :
             other.communications.iterateSchemes(sr)) {
      auto proc = scheme_pair.first;
      auto & other_scheme = scheme_pair.second;
      auto & scheme = communications.createScheme(proc, sr);
      scheme.copy(other_scheme);
    }
  }
}

/* -------------------------------------------------------------------------- */
template <class Entity>
SynchronizerImpl<Entity> & SynchronizerImpl<Entity>::
operator=(const SynchronizerImpl & other) {
  copySchemes(other);
  return *this;
}

/* -------------------------------------------------------------------------- */
template <class Entity>
UInt SynchronizerImpl<Entity>::sanityCheckDataSize(const Array<Entity> &,
                                                   const SynchronizationTag &,
                                                   bool is_comm_desc) const {
  if (not is_comm_desc) {
    return 0;
  }

  UInt size = 0;
  size += sizeof(SynchronizationTag); // tag
  size += sizeof(UInt);               // comm_desc.getNbData();
  size += sizeof(UInt);               // comm_desc.getProc();
  size += sizeof(this->rank);         // mesh.getCommunicator().whoAmI();

  return size;
}

/* -------------------------------------------------------------------------- */
template <class Entity>
void SynchronizerImpl<Entity>::packSanityCheckData(
    CommunicationDescriptor<Entity> & comm_desc) const {
  auto & buffer = comm_desc.getBuffer();
  buffer << comm_desc.getTag();
  buffer << comm_desc.getNbData();
  buffer << comm_desc.getProc();
  buffer << this->rank;

  const auto & tag = comm_desc.getTag();
  const auto & send_element = comm_desc.getScheme();

  this->packSanityCheckData(buffer, send_element, tag);
}

/* -------------------------------------------------------------------------- */
template <class Entity>
void SynchronizerImpl<Entity>::unpackSanityCheckData(
    CommunicationDescriptor<Entity> & comm_desc) const {
  auto & buffer = comm_desc.getBuffer();
  const auto & tag = comm_desc.getTag();

  auto nb_data = comm_desc.getNbData();
  auto proc = comm_desc.getProc();
  auto rank = this->rank;

  decltype(nb_data) recv_nb_data;
  decltype(proc) recv_proc;
  decltype(rank) recv_rank;

  SynchronizationTag t;
  buffer >> t;
  buffer >> recv_nb_data;
  buffer >> recv_proc;
  buffer >> recv_rank;

  AKANTU_DEBUG_ASSERT(
      t == tag, "The tag received does not correspond to the tag expected");

  AKANTU_DEBUG_ASSERT(
      nb_data == recv_nb_data,
      "The nb_data received does not correspond to the nb_data expected");

  AKANTU_DEBUG_ASSERT(UInt(recv_rank) == proc,
                      "The rank received does not correspond to the proc");

  AKANTU_DEBUG_ASSERT(recv_proc == UInt(rank),
                      "The proc received does not correspond to the rank");

  auto & recv_element = comm_desc.getScheme();
  this->unpackSanityCheckData(buffer, recv_element, tag, proc, rank);
}

/* -------------------------------------------------------------------------- */
} // namespace akantu

#endif /* __AKANTU_SYNCHRONIZER_IMPL_TMPL_HH__ */
