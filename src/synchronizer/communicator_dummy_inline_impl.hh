/**
 * @file   communicator_dummy_inline_impl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Nov 07 2017
 * @date last modification: Fri Nov 10 2017
 *
 * @brief  Dummy communicator to make everything work im sequential
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
#include "communicator.hh"
/* -------------------------------------------------------------------------- */
#include <cstring>
#include <type_traits>
#include <vector>
/* -------------------------------------------------------------------------- */

namespace akantu {

Communicator::Communicator(int & /*argc*/, char **& /*argv*/,
                           const private_member & /*unused*/) {}

Communicator::Communicator(const private_member & /*unused*/) {}

template <typename T>
void Communicator::sendImpl(const T * /*unused*/, Int /*unused*/,
                            Int /*unused*/, Int /*unused*/,
                            const CommunicationMode & /*unused*/) const {}
template <typename T>
void Communicator::receiveImpl(T * /*unused*/, Int /*unused*/, Int /*unused*/,
                               Int /*unused*/) const {}

template <typename T>
CommunicationRequest
Communicator::asyncSendImpl(const T * /*unused*/, Int /*unused*/,
                            Int /*unused*/, Int /*unused*/,
                            const CommunicationMode & /*unused*/) const {
  return std::shared_ptr<InternalCommunicationRequest>(
      new InternalCommunicationRequest(0, 0));
}

template <typename T>
CommunicationRequest
Communicator::asyncReceiveImpl(T * /*unused*/, Int /*unused*/, Int /*unused*/,
                               Int /*unused*/) const {
  return std::shared_ptr<InternalCommunicationRequest>(
      new InternalCommunicationRequest(0, 0));
}

template <typename T>
void Communicator::probe(Int /*unused*/, Int /*unused*/,
                         CommunicationStatus & /*unused*/) const {}
template <typename T>
bool Communicator::asyncProbe(Int /*unused*/, Int /*unused*/,
                              CommunicationStatus & /*unused*/) const {
  return true;
}

bool Communicator::test(CommunicationRequest & /*unused*/) { return true; }
bool Communicator::testAll(std::vector<CommunicationRequest> & /*unused*/) {
  return true;
}
void Communicator::wait(CommunicationRequest & /*unused*/) {}
void Communicator::waitAll(std::vector<CommunicationRequest> & /*unused*/) {}
UInt Communicator::waitAny(std::vector<CommunicationRequest> & /*unused*/) {
  return UInt(-1);
}

void Communicator::barrier() const {}
CommunicationRequest Communicator::asyncBarrier() const {
  return std::shared_ptr<InternalCommunicationRequest>(
      new InternalCommunicationRequest(0, 0));
}

template <typename T>
void Communicator::reduceImpl(T * /*unused*/, int /*unused*/,
                              SynchronizerOperation /*unused*/,
                              int /*unused*/) const {}

template <typename T>
void Communicator::allReduceImpl(T * /*unused*/, int /*unused*/,
                                 SynchronizerOperation /*unused*/) const {}

template <typename T>
void Communicator::scanImpl(T * values, T * result, int n,
                            SynchronizerOperation /*unused*/) const {
  if (values == result) {
    return;
  }

  std::copy_n(values, n, result);
}

template <typename T>
void Communicator::exclusiveScanImpl(T * /*values*/, T * result, int n,
                                     SynchronizerOperation /*unused*/) const {
  std::fill_n(result, n, T());
}

template <typename T>
inline void Communicator::allGatherImpl(T * /*unused*/, int /*unused*/) const {}
template <typename T>
inline void Communicator::allGatherVImpl(T * /*unused*/,
                                         int * /*unused*/) const {}

template <typename T>
inline void Communicator::gatherImpl(T * /*unused*/, int /*unused*/,
                                     int /*unused*/) const {}
template <typename T>
void Communicator::gatherImpl(T * values, int nb_values, T * gathered,
                              int /*unused*/) const {
  static_assert(std::is_trivially_copyable<T>{},
                "Cannot send this type of data");
  std::memcpy(gathered, values, nb_values);
}

template <typename T>
inline void Communicator::gatherVImpl(T * /*unused*/, int * /*unused*/,
                                      int /*unused*/) const {}
template <typename T>
inline void Communicator::broadcastImpl(T * /*unused*/, int /*unused*/,
                                        int /*unused*/) const {}

int Communicator::getMaxTag() const { return std::numeric_limits<int>::max(); }
int Communicator::getMinTag() const { return 0; }

Int Communicator::getNbProc() const { return 1; }
Int Communicator::whoAmI() const { return 0; }

} // namespace akantu
