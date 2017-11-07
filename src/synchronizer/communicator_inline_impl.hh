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
#include "communicator.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_STATIC_COMMUNICATOR_INLINE_IMPL_HH__
#define __AKANTU_STATIC_COMMUNICATOR_INLINE_IMPL_HH__

namespace akantu {
/* -------------------------------------------------------------------------- */
inline void
Communicator::freeCommunicationRequest(CommunicationRequest & request) const {
  request.free();
}

/* -------------------------------------------------------------------------- */
inline void Communicator::freeCommunicationRequest(
    std::vector<CommunicationRequest> & requests) const {
  std::vector<CommunicationRequest>::iterator it;
  for (it = requests.begin(); it != requests.end(); ++it) {
    it->free();
  }
}

/* -------------------------------------------------------------------------- */
template <typename T, typename MsgProcessor, typename TagGen>
inline void Communicator::receiveAnyNumber(
    std::vector<CommunicationRequest> & send_requests, Array<T> receive_buffer,
    MsgProcessor && processor, TagGen && tag_gen) const {
  CommunicationRequest barrier_request;
  bool got_all = false, are_send_finished = false;
  while (not got_all) {
    bool are_receives_ready = true;
    while (are_receives_ready) {
      auto && tag = std::forward<TagGen>(tag_gen)();
      CommunicationStatus status;
      are_receives_ready =
          asyncProbe<UInt>(_any_source, tag, status);
      if (are_receives_ready) {
        receive_buffer.resize(status.size());
        receive(receive_buffer, status.getSource(), tag);
        std::forward<MsgProcessor>(processor)(status.getSource(),
                                              receive_buffer);
      }
    }

    if (not are_send_finished) {
      are_send_finished = testAll(send_requests);
      if (are_send_finished)
        barrier_request = asyncBarrier();
    }

    if (are_send_finished) {
      got_all = test(barrier_request);
    }
  }
}
} // namespace akantu

#endif /* __AKANTU_STATIC_COMMUNICATOR_INLINE_IMPL_HH__ */
