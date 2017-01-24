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
} // debug

/* -------------------------------------------------------------------------- */
template <class Entity> class Communications {
public:
  typedef Array<Entity> Scheme;

protected:
  typedef std::map<UInt, Communication> CommunicationPerProcs;
  typedef std::map<SynchronizationTag, CommunicationPerProcs>
      CommunicationsPerTags;
  typedef std::map<UInt, Scheme> CommunicationSchemes;
  typedef std::map<UInt, std::vector<CommunicationRequest> > Request;

  friend class CommunicationDescriptor<Entity>;

public:
  typedef typename CommunicationSchemes::iterator scheme_iterator;
  typedef typename CommunicationSchemes::const_iterator const_scheme_iterator;

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
  void initializeCommunication(const SynchronizationTag & tag, UInt proc, UInt size,
                               const CommunicationSendRecv & sr);
public:
  Communications(const StaticCommunicator & communicator);
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
  UInt incrementCounter(const SynchronizationTag & tag);

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
  const_scheme_iterator begin_send_scheme() const;
  const_scheme_iterator end_send_scheme() const;

  /* ------------------------------------------------------------------------ */
  const_scheme_iterator begin_recv_scheme() const;
  const_scheme_iterator end_recv_scheme() const;

  /* ------------------------------------------------------------------------ */
  scheme_iterator begin_send_scheme();
  scheme_iterator end_send_scheme();

  /* ------------------------------------------------------------------------ */
  scheme_iterator begin_recv_scheme();
  scheme_iterator end_recv_scheme();

  /* ------------------------------------------------------------------------ */
  Scheme & createSendScheme(UInt proc);
  Scheme & createRecvScheme(UInt proc);

  /* ------------------------------------------------------------------------ */
  const Scheme & getSendScheme(UInt proc) const;
  const Scheme & getRecvScheme(UInt proc) const;

  /* ------------------------------------------------------------------------ */
  void resetSchemes();
  /* ------------------------------------------------------------------------ */
  void initializeSendCommunication(const SynchronizationTag & tag, UInt proc,
                                   UInt size);
  void initializeRecvCommunication(const SynchronizationTag & tag, UInt proc,
                                   UInt size);

protected:
  CommunicationSchemes schemes[2];
  CommunicationsPerTags communications[2];
  std::map<SynchronizationTag, UInt> comm_counter;
  std::map<SynchronizationTag, UInt> pending_communications[2];
  const StaticCommunicator & communicator;
};

} // akantu

#include "communications_tmpl.hh"

#endif /* __AKANTU_COMMUNICATIONS_HH__ */
