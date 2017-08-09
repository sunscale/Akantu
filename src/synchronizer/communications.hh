/**
 * @file   communications.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Aug 24 13:56:14 2016
 *
 * @brief
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
#include "communication_descriptor.hh"
/* -------------------------------------------------------------------------- */
#include <map>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_COMMUNICATIONS_HH__
#define __AKANTU_COMMUNICATIONS_HH__

namespace akantu {

namespace debug {
  class CommunicationException : public Exception {
  public:
    CommunicationException()
        : Exception("An exception happen during a communication process.") {}
  };
} // namespace debug

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
  explicit Communications(const StaticCommunicator & communicator);
  /* ------------------------------------------------------------------------ */
  iterator begin_send(const SynchronizationTag & tag);
  iterator end_send(const SynchronizationTag & tag);

  /* ------------------------------------------------------------------------ */
  iterator begin_recv(const SynchronizationTag & tag);
  iterator end_recv(const SynchronizationTag & tag);

  /* ------------------------------------------------------------------------ */
  tag_iterator begin_tag();
  tag_iterator end_tag();

  /* ------------------------------------------------------------------------ */
  bool hasCommunication(const SynchronizationTag & tag) const;
  void incrementCounter(const SynchronizationTag & tag);
  UInt getCounter(const SynchronizationTag & tag) const;

  bool hasPendingRecv(const SynchronizationTag & tag) const;
  bool hasPendingSend(const SynchronizationTag & tag) const;

  const StaticCommunicator & getCommunicator() const;

  /* ------------------------------------------------------------------------ */
  iterator waitAnyRecv(const SynchronizationTag & tag);
  iterator waitAnySend(const SynchronizationTag & tag);

  void waitAllRecv(const SynchronizationTag & tag);
  void waitAllSend(const SynchronizationTag & tag);

  void freeSendRequests(const SynchronizationTag & tag);
  void freeRecvRequests(const SynchronizationTag & tag);

  /* ------------------------------------------------------------------------ */
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
  const StaticCommunicator & communicator;
};

} // namespace akantu

#include "communications_tmpl.hh"

#endif /* __AKANTU_COMMUNICATIONS_HH__ */
