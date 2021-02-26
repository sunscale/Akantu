/**
 * @file   communications.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Sep 07 2016
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Class handling the pending communications and the communications
 * schemes
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
#include "communication_descriptor.hh"
#include "communicator.hh"
/* -------------------------------------------------------------------------- */
#include <map>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_COMMUNICATIONS_HH_
#define AKANTU_COMMUNICATIONS_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class Entity> class Communications {
public:
  using Scheme = Array<Entity>;

protected:
  using CommunicationPerProcs = std::map<UInt, Communication>;
  using CommunicationsPerTags =
      std::map<SynchronizationTag, CommunicationPerProcs>;
  using CommunicationSchemes = std::map<UInt, Scheme>;
  using Request = std::map<UInt, std::vector<CommunicationRequest>>;

  friend class CommunicationDescriptor<Entity>;

public:
  using scheme_iterator = typename CommunicationSchemes::iterator;
  using const_scheme_iterator = typename CommunicationSchemes::const_iterator;

  /* ------------------------------------------------------------------------ */
  class iterator;
  class tag_iterator;
  /* ------------------------------------------------------------------------ */

public:
  CommunicationPerProcs & getCommunications(const SynchronizationTag & tag,
                                            const CommunicationSendRecv & sr);

  /* ------------------------------------------------------------------------ */
  bool hasPending(const SynchronizationTag & tag,
                  const CommunicationSendRecv & sr) const;
  UInt getPending(const SynchronizationTag & tag,
                  const CommunicationSendRecv & sr) const;

  /* ------------------------------------------------------------------------ */
  iterator begin(const SynchronizationTag & tag,
                 const CommunicationSendRecv & sr);
  iterator end(const SynchronizationTag & tag,
               const CommunicationSendRecv & sr);

  /* ------------------------------------------------------------------------ */
  iterator waitAny(const SynchronizationTag & tag,
                   const CommunicationSendRecv & sr);

  /* ------------------------------------------------------------------------ */
  void waitAll(const SynchronizationTag & tag,
               const CommunicationSendRecv & sr);
  void incrementPending(const SynchronizationTag & tag,
                        const CommunicationSendRecv & sr);
  void decrementPending(const SynchronizationTag & tag,
                        const CommunicationSendRecv & sr);
  void freeRequests(const SynchronizationTag & tag,
                    const CommunicationSendRecv & sr);

  /* ------------------------------------------------------------------------ */
  Scheme & createScheme(UInt proc, const CommunicationSendRecv & sr);
  void resetSchemes(const CommunicationSendRecv & sr);

  /* ------------------------------------------------------------------------ */
  void setCommunicationSize(const SynchronizationTag & tag, UInt proc,
                            UInt size, const CommunicationSendRecv & sr);

public:
  explicit Communications(const Communicator & communicator);
  explicit Communications(const Communications & other);

  /* ------------------------------------------------------------------------ */
  void swapSendRecv();

  /* ------------------------------------------------------------------------ */
  class IterableCommunicationDesc {
  public:
    IterableCommunicationDesc(Communications & communications,
                              SynchronizationTag tag, CommunicationSendRecv sr)
        : communications(communications), tag(tag), sr(sr) {}
    auto begin() { return communications.begin(tag, sr); }
    auto end() { return communications.end(tag, sr); }

  private:
    Communications & communications;
    SynchronizationTag tag;
    CommunicationSendRecv sr;
  };

  auto iterateRecv(const SynchronizationTag & tag) {
    return IterableCommunicationDesc(*this, tag, _recv);
  }
  auto iterateSend(const SynchronizationTag & tag) {
    return IterableCommunicationDesc(*this, tag, _send);
  }

  /* ------------------------------------------------------------------------ */
  // iterator begin_send(const SynchronizationTag & tag);
  // iterator end_send(const SynchronizationTag & tag);

  /* ------------------------------------------------------------------------ */
  // iterator begin_recv(const SynchronizationTag & tag);
  // iterator end_recv(const SynchronizationTag & tag);

  /* ------------------------------------------------------------------------ */
  class IterableTags {
  public:
    explicit IterableTags(Communications & communications)
        : communications(communications) {}
    decltype(auto) begin() { return communications.begin_tag(); }
    decltype(auto) end() { return communications.end_tag(); }

  private:
    Communications & communications;
  };

  decltype(auto) iterateTags() { return IterableTags(*this); }

  tag_iterator begin_tag();
  tag_iterator end_tag();

  /* ------------------------------------------------------------------------ */
  bool hasCommunication(const SynchronizationTag & tag) const;
  void incrementCounter(const SynchronizationTag & tag);
  UInt getCounter(const SynchronizationTag & tag) const;
  bool hasCommunicationSize(const SynchronizationTag & tag) const;
  void invalidateSizes();

  /* ------------------------------------------------------------------------ */
  bool hasPendingRecv(const SynchronizationTag & tag) const;
  bool hasPendingSend(const SynchronizationTag & tag) const;

  const auto & getCommunicator() const;

  /* ------------------------------------------------------------------------ */
  iterator waitAnyRecv(const SynchronizationTag & tag);
  iterator waitAnySend(const SynchronizationTag & tag);

  void waitAllRecv(const SynchronizationTag & tag);
  void waitAllSend(const SynchronizationTag & tag);

  void freeSendRequests(const SynchronizationTag & tag);
  void freeRecvRequests(const SynchronizationTag & tag);

  /* ------------------------------------------------------------------------ */
  /* ------------------------------------------------------------------------ */
  class IterableSchemes {
  public:
    IterableSchemes(Communications & communications, CommunicationSendRecv sr)
        : communications(communications), sr(sr) {}
    decltype(auto) begin() { return communications.begin_scheme(sr); }
    decltype(auto) end() { return communications.end_scheme(sr); }

  private:
    Communications & communications;
    CommunicationSendRecv sr;
  };

  class ConstIterableSchemes {
  public:
    ConstIterableSchemes(const Communications & communications,
                         CommunicationSendRecv sr)
        : communications(communications), sr(sr) {}
    decltype(auto) begin() const { return communications.begin_scheme(sr); }
    decltype(auto) end() const { return communications.end_scheme(sr); }

  private:
    const Communications & communications;
    CommunicationSendRecv sr;
  };

  decltype(auto) iterateSchemes(const CommunicationSendRecv & sr) {
    return IterableSchemes(*this, sr);
  }
  decltype(auto) iterateSchemes(const CommunicationSendRecv & sr) const {
    return ConstIterableSchemes(*this, sr);
  }

  decltype(auto) iterateSendSchemes() { return IterableSchemes(*this, _send); }
  decltype(auto) iterateSendSchemes() const {
    return ConstIterableSchemes(*this, _send);
  }

  decltype(auto) iterateRecvSchemes() { return IterableSchemes(*this, _recv); }
  decltype(auto) iterateRecvSchemes() const {
    return ConstIterableSchemes(*this, _recv);
  }

  scheme_iterator begin_scheme(const CommunicationSendRecv & sr);
  scheme_iterator end_scheme(const CommunicationSendRecv & sr);
  const_scheme_iterator begin_scheme(const CommunicationSendRecv & sr) const;
  const_scheme_iterator end_scheme(const CommunicationSendRecv & sr) const;

  /* ------------------------------------------------------------------------ */
  scheme_iterator begin_send_scheme();
  scheme_iterator end_send_scheme();
  const_scheme_iterator begin_send_scheme() const;
  const_scheme_iterator end_send_scheme() const;

  /* ------------------------------------------------------------------------ */
  scheme_iterator begin_recv_scheme();
  scheme_iterator end_recv_scheme();
  const_scheme_iterator begin_recv_scheme() const;
  const_scheme_iterator end_recv_scheme() const;

  /* ------------------------------------------------------------------------ */
  Scheme & createSendScheme(UInt proc);
  Scheme & createRecvScheme(UInt proc);

  /* ------------------------------------------------------------------------ */
  Scheme & getScheme(UInt proc, const CommunicationSendRecv & sr);
  const Scheme & getScheme(UInt proc, const CommunicationSendRecv & sr) const;

  /* ------------------------------------------------------------------------ */
  void resetSchemes();
  /* ------------------------------------------------------------------------ */
  void setSendCommunicationSize(const SynchronizationTag & tag, UInt proc,
                                UInt size);
  void setRecvCommunicationSize(const SynchronizationTag & tag, UInt proc,
                                UInt size);

  void initializeCommunications(const SynchronizationTag & tag);

protected:
  CommunicationSchemes schemes[2];
  CommunicationsPerTags communications[2];

  std::map<SynchronizationTag, UInt> comm_counter;
  std::map<SynchronizationTag, UInt> pending_communications[2];
  std::map<SynchronizationTag, bool> comm_size_computed;

  const Communicator & communicator;
};

} // namespace akantu

#include "communications_tmpl.hh"

#endif /* AKANTU_COMMUNICATIONS_HH_ */
