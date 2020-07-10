/**
 * @file   communications_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Sep 07 2016
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Implementation of Communications
 *
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
#include "communications.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_COMMUNICATIONS_TMPL_HH__
#define __AKANTU_COMMUNICATIONS_TMPL_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class Entity>
Communications<Entity>::Communications(const Communicator & communicator)
    : communicator(communicator) {}

/* -------------------------------------------------------------------------- */
template <class Entity>
Communications<Entity>::Communications(const Communications & other)
    : communicator(other.communicator) {
  for (auto sr : iterate_send_recv) {
    for (const auto & scheme_pair : other.iterateSchemes(sr)) {
      auto proc = scheme_pair.first;
      auto & other_scheme = scheme_pair.second;
      auto & scheme = this->createScheme(proc, sr);
      scheme.copy(other_scheme);
    }
  }

  this->invalidateSizes();
}

/* -------------------------------------------------------------------------- */
template <class Entity> void Communications<Entity>::swapSendRecv() {
  std::swap(schemes[_send], schemes[_recv]);
}

/* -------------------------------------------------------------------------- */
template <class Entity> class Communications<Entity>::iterator {
  using communication_iterator =
      typename std::map<UInt, Communication>::iterator;

public:
  iterator() : communications(nullptr) {}
  iterator(scheme_iterator scheme_it, communication_iterator comm_it,
           Communications<Entity> & communications,
           const SynchronizationTag & tag)
      : scheme_it(scheme_it), comm_it(comm_it), communications(&communications),
        tag(tag) {}

  iterator(const iterator & other) = default;
  iterator(iterator && other) noexcept = default;
  iterator & operator=(const iterator & other) = default;
  iterator & operator=(iterator && other) noexcept = default;

  iterator & operator++() {
    ++scheme_it;
    ++comm_it;
    return *this;
  }

  CommunicationDescriptor<Entity> operator*() {
    AKANTU_DEBUG_ASSERT(
        scheme_it->first == comm_it->first,
        "The two iterators are not in phase, something wrong"
            << " happened, time to take out your favorite debugger ("
            << scheme_it->first << " != " << comm_it->first << ")");
    return CommunicationDescriptor<Entity>(comm_it->second, scheme_it->second,
                                           *communications, tag,
                                           scheme_it->first);
  }

  bool operator==(const iterator & other) const {
    return scheme_it == other.scheme_it && comm_it == other.comm_it;
  }

  bool operator!=(const iterator & other) const {
    return scheme_it != other.scheme_it || comm_it != other.comm_it;
  }

private:
  scheme_iterator scheme_it;
  communication_iterator comm_it;
  Communications<Entity> * communications;
  SynchronizationTag tag;
};

/* -------------------------------------------------------------------------- */
template <class Entity> class Communications<Entity>::tag_iterator {
  using internal_iterator = std::map<SynchronizationTag, UInt>::const_iterator;

public:
  tag_iterator(const internal_iterator & it) : it(it) {}
  tag_iterator & operator++() {
    ++it;
    return *this;
  }
  SynchronizationTag operator*() { return it->first; }
  bool operator==(const tag_iterator & other) const { return it == other.it; }
  bool operator!=(const tag_iterator & other) const { return it != other.it; }

private:
  internal_iterator it;
};

/* -------------------------------------------------------------------------- */
template <class Entity>
typename Communications<Entity>::CommunicationPerProcs &
Communications<Entity>::getCommunications(const SynchronizationTag & tag,
                                          const CommunicationSendRecv & sr) {
  auto comm_it = this->communications[sr].find(tag);
  if (comm_it == this->communications[sr].end())
    AKANTU_CUSTOM_EXCEPTION_INFO(
        debug::CommunicationException(),
        "No known communications for the tag: " << tag);
  return comm_it->second;
}

/* ---------------------------------------------------------------------- */
template <class Entity>
UInt Communications<Entity>::getPending(
    const SynchronizationTag & tag, const CommunicationSendRecv & sr) const {
  const std::map<SynchronizationTag, UInt> & pending =
      pending_communications[sr];
  auto it = pending.find(tag);

  if (it == pending.end())
    return 0;
  return it->second;
}

/* -------------------------------------------------------------------------- */
template <class Entity>
bool Communications<Entity>::hasPending(
    const SynchronizationTag & tag, const CommunicationSendRecv & sr) const {
  return this->hasCommunication(tag) && (this->getPending(tag, sr) != 0);
}

/* ---------------------------------------------------------------------- */
template <class Entity>
typename Communications<Entity>::iterator
Communications<Entity>::begin(const SynchronizationTag & tag,
                              const CommunicationSendRecv & sr) {
  auto & comms = this->getCommunications(tag, sr);
  return iterator(this->schemes[sr].begin(), comms.begin(), *this, tag);
}

template <class Entity>
typename Communications<Entity>::iterator
Communications<Entity>::end(const SynchronizationTag & tag,
                            const CommunicationSendRecv & sr) {
  auto & comms = this->getCommunications(tag, sr);
  return iterator(this->schemes[sr].end(), comms.end(), *this, tag);
}

/* ---------------------------------------------------------------------- */
template <class Entity>
typename Communications<Entity>::iterator
Communications<Entity>::waitAny(const SynchronizationTag & tag,
                                const CommunicationSendRecv & sr) {
  auto & comms = this->getCommunications(tag, sr);
  auto it = comms.begin();
  auto end = comms.end();

  std::vector<CommunicationRequest> requests;

  for (; it != end; ++it) {
    auto & request = it->second.request();
    if (!request.isFreed())
      requests.push_back(request);
  }

  UInt req_id = communicator.waitAny(requests);
  if (req_id != UInt(-1)) {
    auto & request = requests[req_id];
    UInt proc = sr == _recv ? request.getSource() : request.getDestination();

    return iterator(this->schemes[sr].find(proc), comms.find(proc), *this, tag);
  } else {
    return this->end(tag, sr);
  }
}

/* ---------------------------------------------------------------------- */
template <class Entity>
void Communications<Entity>::waitAll(const SynchronizationTag & tag,
                                     const CommunicationSendRecv & sr) {
  auto & comms = this->getCommunications(tag, sr);
  auto it = comms.begin();
  auto end = comms.end();

  std::vector<CommunicationRequest> requests;

  for (; it != end; ++it) {
    requests.push_back(it->second.request());
  }

  communicator.waitAll(requests);
}

template <class Entity>
void Communications<Entity>::incrementPending(
    const SynchronizationTag & tag, const CommunicationSendRecv & sr) {
  ++(pending_communications[sr][tag]);
}

template <class Entity>
void Communications<Entity>::decrementPending(
    const SynchronizationTag & tag, const CommunicationSendRecv & sr) {
  --(pending_communications[sr][tag]);
}

template <class Entity>
void Communications<Entity>::freeRequests(const SynchronizationTag & tag,
                                          const CommunicationSendRecv & sr) {
  iterator it = this->begin(tag, sr);
  iterator end = this->end(tag, sr);

  for (; it != end; ++it) {
    (*it).freeRequest();
  }
}

/* -------------------------------------------------------------------------- */
template <class Entity>
typename Communications<Entity>::Scheme &
Communications<Entity>::createScheme(UInt proc,
                                     const CommunicationSendRecv & sr) {
  // scheme_iterator it = schemes[sr].find(proc);
  // if (it != schemes[sr].end()) {
  //   AKANTU_CUSTOM_EXCEPTION_INFO(debug::CommunicationException(),
  //                                "Communication scheme("
  //                                    << sr
  //                                    << ") already created for proc: " <<
  //                                    proc);
  // }
  return schemes[sr][proc];
}

template <class Entity>
void Communications<Entity>::resetSchemes(const CommunicationSendRecv & sr) {
  auto it = this->schemes[sr].begin();
  auto end = this->schemes[sr].end();
  for (; it != end; ++it) {
    it->second.resize(0);
  }
}

/* -------------------------------------------------------------------------- */
template <class Entity>
void Communications<Entity>::setCommunicationSize(
    const SynchronizationTag & tag, UInt proc, UInt size,
    const CommunicationSendRecv & sr) {
  // accessor that fails if it does not exists
  comm_size_computed[tag] = true; // TODO: need perhaps to be split based on sr
  auto & comms = this->communications[sr];
  auto & comms_per_tag = comms.at(tag);

  comms_per_tag.at(proc).resize(size);
}

/* -------------------------------------------------------------------------- */
template <class Entity>
void Communications<Entity>::initializeCommunications(
    const SynchronizationTag & tag) {

  for (auto t : send_recv_t{}) {
    pending_communications[t].insert(std::make_pair(tag, 0));

    auto & comms = this->communications[t];
    auto & comms_per_tag =
        comms.insert(std::make_pair(tag, CommunicationPerProcs()))
            .first->second;

    for (auto pair : this->schemes[t]) {
      comms_per_tag.emplace(std::piecewise_construct,
                            std::forward_as_tuple(pair.first),
                            std::forward_as_tuple(t));
    }
  }

  comm_counter.insert(std::make_pair(tag, 0));
}

/* -------------------------------------------------------------------------- */
template <class Entity>
typename Communications<Entity>::tag_iterator
Communications<Entity>::begin_tag() {
  return tag_iterator(comm_counter.begin());
}

template <class Entity>
typename Communications<Entity>::tag_iterator
Communications<Entity>::end_tag() {
  return tag_iterator(comm_counter.end());
}

/* -------------------------------------------------------------------------- */
template <class Entity>
typename Communications<Entity>::scheme_iterator
Communications<Entity>::begin_scheme(const CommunicationSendRecv & sr) {
  return this->schemes[sr].begin();
}

template <class Entity>
typename Communications<Entity>::scheme_iterator
Communications<Entity>::end_scheme(const CommunicationSendRecv & sr) {
  return this->schemes[sr].end();
}

/* -------------------------------------------------------------------------- */
template <class Entity>
typename Communications<Entity>::const_scheme_iterator
Communications<Entity>::begin_scheme(const CommunicationSendRecv & sr) const {
  return this->schemes[sr].begin();
}

template <class Entity>
typename Communications<Entity>::const_scheme_iterator
Communications<Entity>::end_scheme(const CommunicationSendRecv & sr) const {
  return this->schemes[sr].end();
}

/* -------------------------------------------------------------------------- */
template <class Entity>
typename Communications<Entity>::scheme_iterator
Communications<Entity>::begin_send_scheme() {
  return this->begin_scheme(_send);
}

template <class Entity>
typename Communications<Entity>::scheme_iterator
Communications<Entity>::end_send_scheme() {
  return this->end_scheme(_send);
}

/* -------------------------------------------------------------------------- */
template <class Entity>
typename Communications<Entity>::const_scheme_iterator
Communications<Entity>::begin_send_scheme() const {
  return this->begin_scheme(_send);
}

template <class Entity>
typename Communications<Entity>::const_scheme_iterator
Communications<Entity>::end_send_scheme() const {
  return this->end_scheme(_send);
}

/* -------------------------------------------------------------------------- */
template <class Entity>
typename Communications<Entity>::scheme_iterator
Communications<Entity>::begin_recv_scheme() {
  return this->begin_scheme(_recv);
}

template <class Entity>
typename Communications<Entity>::scheme_iterator
Communications<Entity>::end_recv_scheme() {
  return this->end_scheme(_recv);
}

/* -------------------------------------------------------------------------- */
template <class Entity>
typename Communications<Entity>::const_scheme_iterator
Communications<Entity>::begin_recv_scheme() const {
  return this->begin_scheme(_recv);
}

template <class Entity>
typename Communications<Entity>::const_scheme_iterator
Communications<Entity>::end_recv_scheme() const {
  return this->end_scheme(_recv);
}

/* ------------------------------------------------------------------------ */
template <class Entity>
bool Communications<Entity>::hasCommunication(
    const SynchronizationTag & tag) const {
  return (communications[_send].find(tag) != communications[_send].end());
}

template <class Entity>
void Communications<Entity>::incrementCounter(const SynchronizationTag & tag) {
  auto it = comm_counter.find(tag);
  if (it == comm_counter.end()) {
    AKANTU_CUSTOM_EXCEPTION_INFO(
        debug::CommunicationException(),
        "No counter initialized in communications for the tags: " << tag);
  }

  ++(it->second);
}

template <class Entity>
UInt Communications<Entity>::getCounter(const SynchronizationTag & tag) const {
  auto it = comm_counter.find(tag);
  if (it == comm_counter.end()) {
    AKANTU_CUSTOM_EXCEPTION_INFO(
        debug::CommunicationException(),
        "No counter initialized in communications for the tags: " << tag);
  }

  return it->second;
}

template <class Entity>
bool Communications<Entity>::hasCommunicationSize(
    const SynchronizationTag & tag) const {
  auto it = comm_size_computed.find(tag);
  if (it == comm_size_computed.end()) {
    return false;
  }

  return it->second;
}

template <class Entity> void Communications<Entity>::invalidateSizes() {
  for (auto && pair : comm_size_computed) {
    pair.second = false;
  }
}

template <class Entity>
bool Communications<Entity>::hasPendingRecv(
    const SynchronizationTag & tag) const {
  return this->hasPending(tag, _recv);
}

template <class Entity>
bool Communications<Entity>::hasPendingSend(
    const SynchronizationTag & tag) const {
  return this->hasPending(tag, _send);
}

template <class Entity>
const auto & Communications<Entity>::getCommunicator() const {
  return communicator;
}

/* -------------------------------------------------------------------------- */
template <class Entity>
typename Communications<Entity>::iterator
Communications<Entity>::waitAnyRecv(const SynchronizationTag & tag) {
  return this->waitAny(tag, _recv);
}

template <class Entity>
typename Communications<Entity>::iterator
Communications<Entity>::waitAnySend(const SynchronizationTag & tag) {
  return this->waitAny(tag, _send);
}

template <class Entity>
void Communications<Entity>::waitAllRecv(const SynchronizationTag & tag) {
  this->waitAll(tag, _recv);
}

template <class Entity>
void Communications<Entity>::waitAllSend(const SynchronizationTag & tag) {
  this->waitAll(tag, _send);
}

template <class Entity>
void Communications<Entity>::freeSendRequests(const SynchronizationTag & tag) {
  this->freeRequests(tag, _send);
}

template <class Entity>
void Communications<Entity>::freeRecvRequests(const SynchronizationTag & tag) {
  this->freeRequests(tag, _recv);
}

/* -------------------------------------------------------------------------- */
template <class Entity>
typename Communications<Entity>::Scheme &
Communications<Entity>::createSendScheme(UInt proc) {
  return createScheme(proc, _send);
}

template <class Entity>
typename Communications<Entity>::Scheme &
Communications<Entity>::createRecvScheme(UInt proc) {
  return createScheme(proc, _recv);
}

/* -------------------------------------------------------------------------- */
template <class Entity> void Communications<Entity>::resetSchemes() {
  resetSchemes(_send);
  resetSchemes(_recv);
}

/* -------------------------------------------------------------------------- */
template <class Entity>
typename Communications<Entity>::Scheme &
Communications<Entity>::getScheme(UInt proc, const CommunicationSendRecv & sr) {
  return this->schemes[sr].find(proc)->second;
}

/* -------------------------------------------------------------------------- */
template <class Entity>
const typename Communications<Entity>::Scheme &
Communications<Entity>::getScheme(UInt proc,
                                  const CommunicationSendRecv & sr) const {
  return this->schemes[sr].find(proc)->second;
}

/* -------------------------------------------------------------------------- */
template <class Entity>
void Communications<Entity>::setSendCommunicationSize(
    const SynchronizationTag & tag, UInt proc, UInt size) {
  this->setCommunicationSize(tag, proc, size, _send);
}

template <class Entity>
void Communications<Entity>::setRecvCommunicationSize(
    const SynchronizationTag & tag, UInt proc, UInt size) {
  this->setCommunicationSize(tag, proc, size, _recv);
}

} // namespace akantu

#endif /* __AKANTU_COMMUNICATIONS_TMPL_HH__ */
