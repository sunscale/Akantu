/**
 * @file   static_communicator_inline_impl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Sep 06 2010
 * @date last modification: Thu Dec 10 2015
 *
 * @brief  implementation of inline functions
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "static_communicator.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_STATIC_COMMUNICATOR_INLINE_IMPL_HH__
#define __AKANTU_STATIC_COMMUNICATOR_INLINE_IMPL_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
inline void
StaticCommunicator::freeCommunicationRequest(CommunicationRequest request) {
  request.free();
}

/* -------------------------------------------------------------------------- */
inline void StaticCommunicator::freeCommunicationRequest(
    std::vector<CommunicationRequest> & requests) {
  std::vector<CommunicationRequest>::iterator it;
  for (it = requests.begin(); it != requests.end(); ++it) {
    it->free();
  }
}

#if defined(__INTEL_COMPILER)
#pragma warning(push)
#pragma warning(disable : 111)
#endif // defined(__INTEL_COMPILER)

/* -------------------------------------------------------------------------- */
#define AKANTU_BOOST_REAL_COMMUNICATOR_CALL(r, call, comm_type)                \
  case BOOST_PP_LIST_AT(comm_type, 0): {                                       \
    BOOST_PP_LIST_AT(comm_type, 1) * comm =                                    \
        static_cast<BOOST_PP_LIST_AT(comm_type, 1) *>(                         \
            real_static_communicator);                                         \
    BOOST_PP_IF(BOOST_PP_LIST_AT(call, 0),                                     \
                return comm->BOOST_PP_LIST_AT(call, 1),                        \
                comm->BOOST_PP_LIST_AT(call, 1);                               \
                break;);                                                       \
  }

#define AKANTU_BOOST_REAL_COMMUNICATOR_SELECT_CALL(call, ret)                  \
  do {                                                                         \
    switch (real_type) {                                                       \
      BOOST_PP_SEQ_FOR_EACH(AKANTU_BOOST_REAL_COMMUNICATOR_CALL,               \
                            (ret, (call, BOST_PP_NIL)),                        \
                            AKANTU_COMMUNICATOR_LIST_ALL)                      \
    default:                                                                   \
      StaticCommunicatorDummy * comm =                                         \
          static_cast<StaticCommunicatorDummy *>(real_static_communicator);    \
      BOOST_PP_IF(ret, return comm->call, comm->call);                         \
    }                                                                          \
  } while (0)

/* -------------------------------------------------------------------------- */
template <typename T>
inline void StaticCommunicator::send(T * buffer, Int size, Int receiver,
                                     Int tag) {
  AKANTU_BOOST_REAL_COMMUNICATOR_SELECT_CALL(send(buffer, size, receiver, tag),
                                             0);
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline void StaticCommunicator::receive(T * buffer, Int size, Int sender,
                                        Int tag) {
  AKANTU_BOOST_REAL_COMMUNICATOR_SELECT_CALL(receive(buffer, size, sender, tag),
                                             0);
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline CommunicationRequest
StaticCommunicator::asyncSend(T * buffer, Int size, Int receiver, Int tag) {
  AKANTU_BOOST_REAL_COMMUNICATOR_SELECT_CALL(
      asyncSend(buffer, size, receiver, tag), 1);
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline CommunicationRequest
StaticCommunicator::asyncReceive(T * buffer, Int size, Int sender, Int tag) {
  AKANTU_BOOST_REAL_COMMUNICATOR_SELECT_CALL(
      asyncReceive(buffer, size, sender, tag), 1);
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline void StaticCommunicator::probe(Int sender, Int tag,
                                      CommunicationStatus & status) {
  AKANTU_BOOST_REAL_COMMUNICATOR_SELECT_CALL(
      template probe<T>(sender, tag, status), 0);
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline void StaticCommunicator::reduce(T * values, int nb_values,
                                       const SynchronizerOperation & op,
                                       int root) {
  AKANTU_BOOST_REAL_COMMUNICATOR_SELECT_CALL(
      reduce(values, nb_values, op, root), 0);
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline void StaticCommunicator::allReduce(T * values, int nb_values,
                                          const SynchronizerOperation & op) {
  AKANTU_BOOST_REAL_COMMUNICATOR_SELECT_CALL(allReduce(values, nb_values, op),
                                             0);
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline void StaticCommunicator::allGather(T * values, int nb_values) {
  AKANTU_BOOST_REAL_COMMUNICATOR_SELECT_CALL(allGather(values, nb_values), 0);
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline void StaticCommunicator::allGatherV(T * values, int * nb_values) {
  AKANTU_BOOST_REAL_COMMUNICATOR_SELECT_CALL(allGatherV(values, nb_values), 0);
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline void StaticCommunicator::gather(T * values, int nb_values, int root) {
  AKANTU_BOOST_REAL_COMMUNICATOR_SELECT_CALL(gather(values, nb_values, root),
                                             0);
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline void StaticCommunicator::gather(T * values, int nb_values, T * gathered,
                                       int nb_gathered) {
  AKANTU_BOOST_REAL_COMMUNICATOR_SELECT_CALL(
      gather(values, nb_values, gathered, nb_gathered), 0);
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline void StaticCommunicator::gatherV(T * values, int * nb_values, int root) {
  AKANTU_BOOST_REAL_COMMUNICATOR_SELECT_CALL(gatherV(values, nb_values, root),
                                             0);
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline void StaticCommunicator::broadcast(T * values, int nb_values, int root) {
  AKANTU_BOOST_REAL_COMMUNICATOR_SELECT_CALL(broadcast(values, nb_values, root),
                                             0);
}

/* -------------------------------------------------------------------------- */
inline void StaticCommunicator::barrier() {
  AKANTU_BOOST_REAL_COMMUNICATOR_SELECT_CALL(barrier(), 0);
}

/* -------------------------------------------------------------------------- */
inline bool StaticCommunicator::testRequest(CommunicationRequest request) {
  AKANTU_BOOST_REAL_COMMUNICATOR_SELECT_CALL(testRequest(request), 1);
}

/* -------------------------------------------------------------------------- */
inline void StaticCommunicator::wait(CommunicationRequest request) {
  AKANTU_BOOST_REAL_COMMUNICATOR_SELECT_CALL(wait(request), 0);
}

/* -------------------------------------------------------------------------- */
inline void
StaticCommunicator::waitAll(std::vector<CommunicationRequest> & requests) {
  AKANTU_BOOST_REAL_COMMUNICATOR_SELECT_CALL(waitAll(requests), 0);
}

/* -------------------------------------------------------------------------- */
inline UInt
StaticCommunicator::waitAny(std::vector<CommunicationRequest> & requests) {
  AKANTU_BOOST_REAL_COMMUNICATOR_SELECT_CALL(waitAny(requests), 1);
}

/* -------------------------------------------------------------------------- */
inline int StaticCommunicator::getMaxTag() {
  AKANTU_BOOST_REAL_COMMUNICATOR_SELECT_CALL(getMaxTag(), 1);
}

/* -------------------------------------------------------------------------- */
inline int StaticCommunicator::getMinTag() {
  AKANTU_BOOST_REAL_COMMUNICATOR_SELECT_CALL(getMinTag(), 1);
}

#if defined(__INTEL_COMPILER)
#pragma warning(pop)
#endif // defined(__INTEL_COMPILER)

}  // akantu

#endif /* __AKANTU_STATIC_COMMUNICATOR_INLINE_IMPL_HH__ */
