/**
 * @file   static_communicator_mpi.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sun Sep 26 2010
 * @date last modification: Thu Jan 21 2016
 *
 * @brief  StaticCommunicatorMPI implementation
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
#include "aka_iterators.hh"
#include "communicator.hh"
#include "mpi_communicator_data.hh"
/* -------------------------------------------------------------------------- */
#include <type_traits>
#include <unordered_map>
#include <vector>
#include <memory>
/* -------------------------------------------------------------------------- */
#include <mpi.h>
/* -------------------------------------------------------------------------- */

#if (defined(__GNUC__) || defined(__GNUG__))
#define GCC_VERSION                                                            \
  (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#if GCC_VERSION < 60000
namespace std {
template <> struct hash<akantu::SynchronizerOperation> {
  using argument_type = akantu::SynchronizerOperation;
  size_t operator()(const argument_type & e) const noexcept {
    auto ue = underlying_type_t<argument_type>(e);
      return uh(ue);
    }

  private:
    const hash<underlying_type_t<argument_type>> uh{};
  };
} // namespace std
#endif
#endif


namespace akantu {

class CommunicationRequestMPI : public InternalCommunicationRequest {
public:
  CommunicationRequestMPI(UInt source, UInt dest)
      : InternalCommunicationRequest(source, dest),
        request(std::make_unique<MPI_Request>()) {}
  MPI_Request & getMPIRequest() { return *request; };

private:
  std::unique_ptr<MPI_Request> request;
};

namespace {
  template <typename T> inline MPI_Datatype getMPIDatatype();
  MPI_Op getMPISynchronizerOperation(const SynchronizerOperation & op) {
    std::unordered_map<SynchronizerOperation, MPI_Op> _operations{
        {SynchronizerOperation::_sum, MPI_SUM},
        {SynchronizerOperation::_min, MPI_MIN},
        {SynchronizerOperation::_max, MPI_MAX},
        {SynchronizerOperation::_prod, MPI_PROD},
        {SynchronizerOperation::_land, MPI_LAND},
        {SynchronizerOperation::_band, MPI_BAND},
        {SynchronizerOperation::_lor, MPI_LOR},
        {SynchronizerOperation::_bor, MPI_BOR},
        {SynchronizerOperation::_lxor, MPI_LXOR},
        {SynchronizerOperation::_bxor, MPI_BXOR},
        {SynchronizerOperation::_min_loc, MPI_MINLOC},
        {SynchronizerOperation::_max_loc, MPI_MAXLOC},
        {SynchronizerOperation::_null, MPI_OP_NULL}};
    return _operations[op];
  }

  template <typename T> MPI_Datatype inline getMPIDatatype() {
    return MPI_DATATYPE_NULL;
  }

#define SPECIALIZE_MPI_DATATYPE(type, mpi_type)                                \
  template <> MPI_Datatype inline getMPIDatatype<type>() { return mpi_type; }

#define COMMA ,
  SPECIALIZE_MPI_DATATYPE(char, MPI_CHAR)
  SPECIALIZE_MPI_DATATYPE(float, MPI_FLOAT)
  SPECIALIZE_MPI_DATATYPE(double, MPI_DOUBLE)
  SPECIALIZE_MPI_DATATYPE(long double, MPI_LONG_DOUBLE)
  SPECIALIZE_MPI_DATATYPE(signed int, MPI_INT)
  SPECIALIZE_MPI_DATATYPE(NodeType, getMPIDatatype<Int>())
  SPECIALIZE_MPI_DATATYPE(unsigned int, MPI_UNSIGNED)
  SPECIALIZE_MPI_DATATYPE(signed long int, MPI_LONG)
  SPECIALIZE_MPI_DATATYPE(unsigned long int, MPI_UNSIGNED_LONG)
  SPECIALIZE_MPI_DATATYPE(signed long long int, MPI_LONG_LONG)
  SPECIALIZE_MPI_DATATYPE(unsigned long long int, MPI_UNSIGNED_LONG_LONG)
  SPECIALIZE_MPI_DATATYPE(SCMinMaxLoc<double COMMA int>, MPI_DOUBLE_INT)
  SPECIALIZE_MPI_DATATYPE(SCMinMaxLoc<float COMMA int>, MPI_FLOAT_INT)
  SPECIALIZE_MPI_DATATYPE(bool, MPI_CXX_BOOL)

  inline int getMPISource(int src) {
    if (src == _any_source)
      return MPI_ANY_SOURCE;
    return src;
  }

  decltype(auto) convertRequests(std::vector<CommunicationRequest> & requests) {
    std::vector<MPI_Request> mpi_requests(requests.size());

    for (auto && request_pair : zip(requests, mpi_requests)) {
      auto && req = std::get<0>(request_pair);
      auto && mpi_req = std::get<1>(request_pair);
      mpi_req = dynamic_cast<CommunicationRequestMPI &>(req.getInternal())
                    .getMPIRequest();
    }
    return mpi_requests;
  }

} // namespace

// this is ugly but shorten the code a lot
#define MPIDATA                                                                \
  (*reinterpret_cast<MPICommunicatorData *>(communicator_data.get()))

/* -------------------------------------------------------------------------- */
/* Implementation                                                             */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
Communicator::Communicator(int & /*argc*/, char **& /*argv*/,
                           const private_member & /*unused*/)
    : communicator_data(std::make_unique<MPICommunicatorData>()) {
  prank = MPIDATA.rank();
  psize = MPIDATA.size();
}

/* -------------------------------------------------------------------------- */
template <typename T>
void Communicator::sendImpl(const T * buffer, Int size, Int receiver, Int tag,
                            const CommunicationMode & mode) const {
  MPI_Comm communicator = MPIDATA.getMPICommunicator();
  MPI_Datatype type = getMPIDatatype<T>();

  switch (mode) {
  case CommunicationMode::_auto:
    MPI_Send(buffer, size, type, receiver, tag, communicator);
    break;
  case CommunicationMode::_synchronous:
    MPI_Ssend(buffer, size, type, receiver, tag, communicator);
    break;
  case CommunicationMode::_ready:
    MPI_Rsend(buffer, size, type, receiver, tag, communicator);
    break;
  }
}

/* -------------------------------------------------------------------------- */
template <typename T>
void Communicator::receiveImpl(T * buffer, Int size, Int sender,
                               Int tag) const {
  MPI_Comm communicator = MPIDATA.getMPICommunicator();
  MPI_Status status;
  MPI_Datatype type = getMPIDatatype<T>();
  MPI_Recv(buffer, size, type, getMPISource(sender), tag, communicator,
           &status);
}

/* -------------------------------------------------------------------------- */
template <typename T>
CommunicationRequest
Communicator::asyncSendImpl(const T * buffer, Int size, Int receiver, Int tag,
                            const CommunicationMode & mode) const {
  MPI_Comm communicator = MPIDATA.getMPICommunicator();
  CommunicationRequestMPI * request =
      new CommunicationRequestMPI(prank, receiver);
  MPI_Request & req = request->getMPIRequest();

  MPI_Datatype type = getMPIDatatype<T>();

  switch (mode) {
  case CommunicationMode::_auto:
    MPI_Isend(buffer, size, type, receiver, tag, communicator, &req);
    break;
  case CommunicationMode::_synchronous:
    MPI_Issend(buffer, size, type, receiver, tag, communicator, &req);
    break;
  case CommunicationMode::_ready:
    MPI_Irsend(buffer, size, type, receiver, tag, communicator, &req);
    break;
  }
  return std::shared_ptr<InternalCommunicationRequest>(request);
}

/* -------------------------------------------------------------------------- */
template <typename T>
CommunicationRequest Communicator::asyncReceiveImpl(T * buffer, Int size,
                                                    Int sender, Int tag) const {
  MPI_Comm communicator = MPIDATA.getMPICommunicator();
  CommunicationRequestMPI * request =
      new CommunicationRequestMPI(sender, prank);
  MPI_Datatype type = getMPIDatatype<T>();

  MPI_Request & req = request->getMPIRequest();
  MPI_Irecv(buffer, size, type, getMPISource(sender), tag, communicator, &req);
  return std::shared_ptr<InternalCommunicationRequest>(request);
}

/* -------------------------------------------------------------------------- */
template <typename T>
void Communicator::probe(Int sender, Int tag,
                         CommunicationStatus & status) const {
  MPI_Comm communicator = MPIDATA.getMPICommunicator();
  MPI_Status mpi_status;
  MPI_Probe(getMPISource(sender), tag, communicator, &mpi_status);

  MPI_Datatype type = getMPIDatatype<T>();
  int count;
  MPI_Get_count(&mpi_status, type, &count);

  status.setSource(mpi_status.MPI_SOURCE);
  status.setTag(mpi_status.MPI_TAG);
  status.setSize(count);
}

/* -------------------------------------------------------------------------- */
template <typename T>
bool Communicator::asyncProbe(Int sender, Int tag,
                              CommunicationStatus & status) const {
  MPI_Comm communicator = MPIDATA.getMPICommunicator();
  MPI_Status mpi_status;
  int test;
  MPI_Iprobe(getMPISource(sender), tag, communicator, &test, &mpi_status);

  if (not test)
    return false;

  MPI_Datatype type = getMPIDatatype<T>();
  int count;
  MPI_Get_count(&mpi_status, type, &count);

  status.setSource(mpi_status.MPI_SOURCE);
  status.setTag(mpi_status.MPI_TAG);
  status.setSize(count);
  return true;
}

/* -------------------------------------------------------------------------- */
bool Communicator::test(CommunicationRequest & request) const {
  MPI_Status status;
  int flag;
  CommunicationRequestMPI & req_mpi =
      dynamic_cast<CommunicationRequestMPI &>(request.getInternal());
  MPI_Request & req = req_mpi.getMPIRequest();

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
      MPI_Test(&req, &flag, &status);

  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Test.");
  return (flag != 0);
}

/* -------------------------------------------------------------------------- */
bool Communicator::testAll(std::vector<CommunicationRequest> & requests) const {
  int are_finished;
  auto && mpi_requests = convertRequests(requests);
  MPI_Testall(mpi_requests.size(), mpi_requests.data(), &are_finished,
              MPI_STATUSES_IGNORE);
  return are_finished != 0 ? true : false;
}

/* -------------------------------------------------------------------------- */
void Communicator::wait(CommunicationRequest & request) const {
  MPI_Status status;
  CommunicationRequestMPI & req_mpi =
      dynamic_cast<CommunicationRequestMPI &>(request.getInternal());
  MPI_Request & req = req_mpi.getMPIRequest();
  MPI_Wait(&req, &status);
}

/* -------------------------------------------------------------------------- */
void Communicator::waitAll(std::vector<CommunicationRequest> & requests) const {
  auto && mpi_requests = convertRequests(requests);
  MPI_Waitall(mpi_requests.size(), mpi_requests.data(), MPI_STATUSES_IGNORE);
}

/* -------------------------------------------------------------------------- */
UInt Communicator::waitAny(std::vector<CommunicationRequest> & requests) const {
  auto && mpi_requests = convertRequests(requests);

  int pos;
  MPI_Waitany(mpi_requests.size(), mpi_requests.data(), &pos,
              MPI_STATUSES_IGNORE);

  if (pos != MPI_UNDEFINED) {
    return pos;
  } else {
    return UInt(-1);
  }
}

/* -------------------------------------------------------------------------- */
void Communicator::barrier() const {
  MPI_Comm communicator = MPIDATA.getMPICommunicator();
  MPI_Barrier(communicator);
}

/* -------------------------------------------------------------------------- */
CommunicationRequest Communicator::asyncBarrier() const {
  MPI_Comm communicator = MPIDATA.getMPICommunicator();
  CommunicationRequestMPI * request = new CommunicationRequestMPI(0, 0);

  MPI_Request & req = request->getMPIRequest();
  MPI_Ibarrier(communicator, &req);

  return std::shared_ptr<InternalCommunicationRequest>(request);
}

/* -------------------------------------------------------------------------- */
template <typename T>
void Communicator::reduceImpl(T * values, int nb_values,
                              const SynchronizerOperation & op,
                              int root) const {
  MPI_Comm communicator = MPIDATA.getMPICommunicator();
  MPI_Datatype type = getMPIDatatype<T>();

  MPI_Reduce(MPI_IN_PLACE, values, nb_values, type,
             getMPISynchronizerOperation(op), root, communicator);
}

/* -------------------------------------------------------------------------- */
template <typename T>
void Communicator::allReduceImpl(T * values, int nb_values,
                                 const SynchronizerOperation & op) const {
  MPI_Comm communicator = MPIDATA.getMPICommunicator();
  MPI_Datatype type = getMPIDatatype<T>();

  MPI_Allreduce(MPI_IN_PLACE, values, nb_values, type,
                getMPISynchronizerOperation(op), communicator);
}

/* -------------------------------------------------------------------------- */
template <typename T>
void Communicator::allGatherImpl(T * values, int nb_values) const {
  MPI_Comm communicator = MPIDATA.getMPICommunicator();
  MPI_Datatype type = getMPIDatatype<T>();

  MPI_Allgather(MPI_IN_PLACE, nb_values, type, values, nb_values, type,
                communicator);
}

/* -------------------------------------------------------------------------- */
template <typename T>
void Communicator::allGatherVImpl(T * values, int * nb_values) const {
  MPI_Comm communicator = MPIDATA.getMPICommunicator();
  std::vector<int> displs(psize);
  displs[0] = 0;
  for (int i = 1; i < psize; ++i) {
    displs[i] = displs[i - 1] + nb_values[i - 1];
  }

  MPI_Datatype type = getMPIDatatype<T>();
  MPI_Allgatherv(MPI_IN_PLACE, *nb_values, type, values, nb_values,
                 displs.data(), type, communicator);
}

/* -------------------------------------------------------------------------- */
template <typename T>
void Communicator::gatherImpl(T * values, int nb_values, int root) const {
  MPI_Comm communicator = MPIDATA.getMPICommunicator();
  T *send_buf = nullptr, *recv_buf = nullptr;
  if (prank == root) {
    send_buf = (T *)MPI_IN_PLACE;
    recv_buf = values;
  } else {
    send_buf = values;
  }

  MPI_Datatype type = getMPIDatatype<T>();
  MPI_Gather(send_buf, nb_values, type, recv_buf, nb_values, type, root,
             communicator);
}

/* -------------------------------------------------------------------------- */
template <typename T>
void Communicator::gatherImpl(T * values, int nb_values, T * gathered,
                              int nb_gathered) const {
  MPI_Comm communicator = MPIDATA.getMPICommunicator();
  T * send_buf = values;
  T * recv_buf = gathered;

  if (nb_gathered == 0)
    nb_gathered = nb_values;

  MPI_Datatype type = getMPIDatatype<T>();
  MPI_Gather(send_buf, nb_values, type, recv_buf, nb_gathered, type,
             this->prank, communicator);
}

/* -------------------------------------------------------------------------- */
template <typename T>
void Communicator::gatherVImpl(T * values, int * nb_values, int root) const {
  MPI_Comm communicator = MPIDATA.getMPICommunicator();
  int * displs = nullptr;
  if (prank == root) {
    displs = new int[psize];
    displs[0] = 0;
    for (int i = 1; i < psize; ++i) {
      displs[i] = displs[i - 1] + nb_values[i - 1];
    }
  }

  T *send_buf = nullptr, *recv_buf = nullptr;
  if (prank == root) {
    send_buf = (T *)MPI_IN_PLACE;
    recv_buf = values;
  } else
    send_buf = values;

  MPI_Datatype type = getMPIDatatype<T>();

  MPI_Gatherv(send_buf, *nb_values, type, recv_buf, nb_values, displs, type,
              root, communicator);

  if (prank == root) {
    delete[] displs;
  }
}

/* -------------------------------------------------------------------------- */
template <typename T>
void Communicator::broadcastImpl(T * values, int nb_values, int root) const {
  MPI_Comm communicator = MPIDATA.getMPICommunicator();
  MPI_Datatype type = getMPIDatatype<T>();
  MPI_Bcast(values, nb_values, type, root, communicator);
}

/* -------------------------------------------------------------------------- */
int Communicator::getMaxTag() const { return MPIDATA.getMaxTag(); }
int Communicator::getMinTag() const { return 0; }

/* -------------------------------------------------------------------------- */

} // namespace akantu
