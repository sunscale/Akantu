/**
 * @file   communication_descriptor_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Sep 07 2016
 * @date last modification: Thu Jan 25 2018
 *
 * @brief  implementation of CommunicationDescriptor
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
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_COMMUNICATION_DESCRIPTOR_TMPL_HH_
#define AKANTU_COMMUNICATION_DESCRIPTOR_TMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
/* Implementations                                                            */
/* -------------------------------------------------------------------------- */
template <class Entity>
CommunicationDescriptor<Entity>::CommunicationDescriptor(
    Communication & communication, Array<Entity> & scheme,
    Communications<Entity> & communications, const SynchronizationTag & tag,
    UInt proc)
    : communication(communication), scheme(scheme),
      communications(communications), tag(tag), proc(proc),
      rank(communications.getCommunicator().whoAmI()) {
  counter = communications.getCounter(tag);
}

/* -------------------------------------------------------------------------- */
template <class Entity>
CommunicationBuffer & CommunicationDescriptor<Entity>::getBuffer() {
  return communication.buffer();
}

/* -------------------------------------------------------------------------- */
template <class Entity>
CommunicationRequest & CommunicationDescriptor<Entity>::getRequest() {
  return communication.request();
}

/* -------------------------------------------------------------------------- */
template <class Entity> void CommunicationDescriptor<Entity>::freeRequest() {
  const auto & comm = communications.getCommunicator();

  // comm.test(communication.request());
  comm.freeCommunicationRequest(communication.request());

  communications.decrementPending(tag, communication.type());
}

/* -------------------------------------------------------------------------- */
template <class Entity>
const Array<Entity> & CommunicationDescriptor<Entity>::getScheme() {
  return scheme;
}

template <class Entity> void CommunicationDescriptor<Entity>::resetBuffer() {
  this->communication.buffer().reset();
}

/* -------------------------------------------------------------------------- */
template <class Entity>
void CommunicationDescriptor<Entity>::packData(
    const DataAccessor<Entity> & accessor) {
  AKANTU_DEBUG_ASSERT(
      communication.type() == _send,
      "Cannot pack data on communication that is not of type _send");
  accessor.packData(communication.buffer(), scheme, tag);
}

/* -------------------------------------------------------------------------- */
template <class Entity>
void CommunicationDescriptor<Entity>::unpackData(
    DataAccessor<Entity> & accessor) {
  AKANTU_DEBUG_ASSERT(
      communication.type() == _recv,
      "Cannot unpack data from communication that is not of type _recv");
  accessor.unpackData(communication.buffer(), scheme, tag);
}

/* -------------------------------------------------------------------------- */
template <class Entity>
void CommunicationDescriptor<Entity>::postSend(int hash_id) {
  AKANTU_DEBUG_ASSERT(communication.type() == _send,
                      "Cannot send a communication that is not of type _send");
  Tag comm_tag = Tag::genTag(rank, counter, tag, hash_id);
  AKANTU_DEBUG_ASSERT(communication.buffer().getPackedSize() ==
                          communication.size(),
                      "a problem have been introduced with "
                          << "false sent sizes declaration "
                          << communication.buffer().getPackedSize()
                          << " != " << communication.size());

  AKANTU_DEBUG_INFO("Posting send to proc " << proc << " (tag: " << tag << " - "
                                            << communication.size()
                                            << " data to send) "
                                            << " [ " << comm_tag << " ]");

  CommunicationRequest & request = communication.request();
  request = communications.getCommunicator().asyncSend(communication.buffer(),
                                                       proc, comm_tag);
  communications.incrementPending(tag, communication.type());
}

/* -------------------------------------------------------------------------- */
template <class Entity>
void CommunicationDescriptor<Entity>::postRecv(int hash_id) {
  AKANTU_DEBUG_ASSERT(communication.type() == _recv,
                      "Cannot receive data for communication ("
                          << communication.type()
                          << ")that is not of type _recv");

  Tag comm_tag = Tag::genTag(proc, counter, tag, hash_id);
  AKANTU_DEBUG_INFO("Posting receive from proc "
                    << proc << " (tag: " << tag << " - " << communication.size()
                    << " data to receive) "
                    << " [ " << comm_tag << " ]");

  CommunicationRequest & request = communication.request();
  request = communications.getCommunicator().asyncReceive(
      communication.buffer(), proc, comm_tag);
  communications.incrementPending(tag, communication.type());
}

} // namespace akantu

#endif /* AKANTU_COMMUNICATION_DESCRIPTOR_TMPL_HH_ */
