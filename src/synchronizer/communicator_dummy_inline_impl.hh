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
void Communicator::sendImpl(const T *, Int, Int, Int,
                            const CommunicationMode &) const {}
template <typename T>
void Communicator::receiveImpl(T *, Int, Int, Int) const {}

template <typename T>
CommunicationRequest
Communicator::asyncSendImpl(const T *, Int, Int, Int,
                            const CommunicationMode &) const {
  return std::shared_ptr<InternalCommunicationRequest>(
      new InternalCommunicationRequest(0, 0));
}

template <typename T>
CommunicationRequest Communicator::asyncReceiveImpl(T *, Int, Int, Int) const {
  return std::shared_ptr<InternalCommunicationRequest>(
      new InternalCommunicationRequest(0, 0));
}

template <typename T>
void Communicator::probe(Int, Int, CommunicationStatus &) const {}
template <typename T>
bool Communicator::asyncProbe(Int, Int, CommunicationStatus &) const {
  return true;
}

bool Communicator::test(CommunicationRequest &) const { return true; }
bool Communicator::testAll(std::vector<CommunicationRequest> &) const {
  return true;
}
void Communicator::wait(CommunicationRequest &) const {}
void Communicator::waitAll(std::vector<CommunicationRequest> &) const {}
UInt Communicator::waitAny(std::vector<CommunicationRequest> &) const {
  return UInt(-1);
}

void Communicator::barrier() const {}
CommunicationRequest Communicator::asyncBarrier() const {
  return std::shared_ptr<InternalCommunicationRequest>(
      new InternalCommunicationRequest(0, 0));
}

template <typename T>
void Communicator::reduceImpl(T *, int, SynchronizerOperation, int) const {}

template <typename T>
void Communicator::allReduceImpl(T *, int, SynchronizerOperation) const {}

template <typename T>
void Communicator::scanImpl(T * values, T * result, int n,
                            SynchronizerOperation) const {
  if (values == result)
    return;

  std::copy_n(values, n, result);
}

template <typename T>
void Communicator::exclusiveScanImpl(T * /*values*/, T * result, int n,
                                     SynchronizerOperation) const {
  std::fill_n(result, n, T());
}

template <typename T> inline void Communicator::allGatherImpl(T *, int) const {}
template <typename T>
inline void Communicator::allGatherVImpl(T *, int *) const {}

template <typename T>
inline void Communicator::gatherImpl(T *, int, int) const {}
template <typename T>
void Communicator::gatherImpl(T * values, int nb_values, T * gathered,
                              int) const {
  static_assert(std::is_trivially_copyable<T>{},
                "Cannot send this type of data");
  std::memcpy(gathered, values, nb_values);
}

template <typename T>
inline void Communicator::gatherVImpl(T *, int *, int) const {}
template <typename T>
inline void Communicator::broadcastImpl(T *, int, int) const {}

int Communicator::getMaxTag() const { return std::numeric_limits<int>::max(); }
int Communicator::getMinTag() const { return 0; }

Int Communicator::getNbProc() const { return 1; }
Int Communicator::whoAmI() const { return 0; }

} // namespace akantu
