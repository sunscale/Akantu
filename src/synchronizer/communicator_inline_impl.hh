/**
 * @file   communicator_inline_impl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Feb 02 2016
 * @date last modification: Tue Nov 07 2017
 *
 * @brief  implementation of inline functions
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
template <typename T, typename MsgProcessor>
inline void Communicator::receiveAnyNumber(
    std::vector<CommunicationRequest> & send_requests, Array<T> receive_buffer,
    MsgProcessor && processor, Int tag) const {
  CommunicationRequest barrier_request;
  bool got_all = false, are_send_finished = false;
  while (not got_all) {
    bool are_receives_ready = true;
    while (are_receives_ready) {
      CommunicationStatus status;
      are_receives_ready = asyncProbe<UInt>(_any_source, tag, status);
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
