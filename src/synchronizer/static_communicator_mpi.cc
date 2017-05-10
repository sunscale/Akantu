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
#include <iostream>

/* -------------------------------------------------------------------------- */
#include "mpi_type_wrapper.hh"
#include "static_communicator_mpi.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

MPI_Op MPITypeWrapper::synchronizer_operation_to_mpi_op[_so_null + 1] = {
    MPI_SUM, MPI_MIN,  MPI_MAX,  MPI_PROD,   MPI_LAND,   MPI_BAND,   MPI_LOR,
    MPI_BOR, MPI_LXOR, MPI_BXOR, MPI_MINLOC, MPI_MAXLOC, MPI_OP_NULL};

class CommunicationRequestMPI : public InternalCommunicationRequest {
public:
  CommunicationRequestMPI(UInt source, UInt dest);
  MPI_Request & getMPIRequest() { return *request; };

private:
  std::unique_ptr<MPI_Request> request;
};

/* -------------------------------------------------------------------------- */
/* Implementation                                                             */
/* -------------------------------------------------------------------------- */
CommunicationRequestMPI::CommunicationRequestMPI(UInt source, UInt dest)
    : InternalCommunicationRequest(source, dest), request(new MPI_Request) {}

/* -------------------------------------------------------------------------- */
StaticCommunicatorMPI::StaticCommunicatorMPI(int & argc, char **& argv)
    : RealStaticCommunicator(argc, argv) {

  int is_initialized = false;
  MPI_Initialized(&is_initialized);
  if (!is_initialized) {
    MPI_Init(&argc, &argv);
  }

  this->is_externaly_initialized = is_initialized;

  mpi_data = new MPITypeWrapper(*this);
  mpi_data->setMPICommunicator(MPI_COMM_WORLD);
}

/* -------------------------------------------------------------------------- */
StaticCommunicatorMPI::~StaticCommunicatorMPI() {
  if (!this->is_externaly_initialized) {
    MPI_Finalize();
    delete this->mpi_data;
  }
}

/* -------------------------------------------------------------------------- */
template <typename T>
void StaticCommunicatorMPI::send(T * buffer, Int size, Int receiver, Int tag) {
  MPI_Comm communicator = mpi_data->getMPICommunicator();
  MPI_Datatype type = MPITypeWrapper::getMPIDatatype<T>();
#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
      MPI_Send(buffer, size, type, receiver, tag, communicator);

  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Send.");
}

/* -------------------------------------------------------------------------- */
template <typename T>
void StaticCommunicatorMPI::receive(T * buffer, Int size, Int sender, Int tag) {
  MPI_Comm communicator = mpi_data->getMPICommunicator();
  MPI_Status status;
  MPI_Datatype type = MPITypeWrapper::getMPIDatatype<T>();

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
      MPI_Recv(buffer, size, type, sender, tag, communicator, &status);

  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Recv.");
}

/* -------------------------------------------------------------------------- */
template <typename T>
CommunicationRequest StaticCommunicatorMPI::asyncSend(T * buffer, Int size,
                                                      Int receiver, Int tag) {
  MPI_Comm communicator = mpi_data->getMPICommunicator();
  CommunicationRequestMPI * request =
      new CommunicationRequestMPI(prank, receiver);
  MPI_Datatype type = MPITypeWrapper::getMPIDatatype<T>();

  MPI_Request & req = request->getMPIRequest();
#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
      MPI_Isend(buffer, size, type, receiver, tag, communicator, &req);

  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Isend.");
  return std::shared_ptr<InternalCommunicationRequest>(request);
}

/* -------------------------------------------------------------------------- */
template <typename T>
CommunicationRequest StaticCommunicatorMPI::asyncReceive(T * buffer, Int size,
                                                         Int sender, Int tag) {
  MPI_Comm communicator = mpi_data->getMPICommunicator();
  CommunicationRequestMPI * request =
      new CommunicationRequestMPI(sender, prank);
  MPI_Datatype type = MPITypeWrapper::getMPIDatatype<T>();

  MPI_Request & req = request->getMPIRequest();
#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
      MPI_Irecv(buffer, size, type, sender, tag, communicator, &req);

  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Irecv.");
  return std::shared_ptr<InternalCommunicationRequest>(request);
}

/* -------------------------------------------------------------------------- */
template <typename T>
void StaticCommunicatorMPI::probe(Int sender, Int tag,
                                  CommunicationStatus & status) {
  MPI_Comm communicator = mpi_data->getMPICommunicator();
  MPI_Status mpi_status;
#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
      MPI_Probe(sender, tag, communicator, &mpi_status);
  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Probe.");

  MPI_Datatype type = MPITypeWrapper::getMPIDatatype<T>();
  int count;
  MPI_Get_count(&mpi_status, type, &count);

  status.setSource(mpi_status.MPI_SOURCE);
  status.setTag(mpi_status.MPI_TAG);
  status.setSize(count);
}

/* -------------------------------------------------------------------------- */
bool StaticCommunicatorMPI::testRequest(CommunicationRequest & request) {
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
void StaticCommunicatorMPI::wait(CommunicationRequest & request) {
  MPI_Status status;
  CommunicationRequestMPI & req_mpi =
      dynamic_cast<CommunicationRequestMPI &>(request.getInternal());
  MPI_Request & req = req_mpi.getMPIRequest();

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
      MPI_Wait(&req, &status);

  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Wait.");
}

/* -------------------------------------------------------------------------- */
void StaticCommunicatorMPI::waitAll(
    std::vector<CommunicationRequest> & requests) {
  MPI_Status status;
  std::vector<CommunicationRequest>::iterator it;
  for (it = requests.begin(); it != requests.end(); ++it) {
    MPI_Request & req =
        dynamic_cast<CommunicationRequestMPI &>(it->getInternal())
            .getMPIRequest();

#if !defined(AKANTU_NDEBUG)
    int ret =
#endif
        MPI_Wait(&req, &status);

    AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Wait.");
  }
}

/* -------------------------------------------------------------------------- */
UInt StaticCommunicatorMPI::waitAny(
    std::vector<CommunicationRequest> & requests) {
  MPI_Status status;
  std::vector<MPI_Request> reqs(requests.size());
  UInt r = 0;
  for (auto it = requests.begin(); it != requests.end(); ++it, ++r) {
    reqs[r] = static_cast<CommunicationRequestMPI *>(&it->getInternal())
                  ->getMPIRequest();
  }

  int pos;
#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
      MPI_Waitany(requests.size(), reqs.data(), &pos, &status);
  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Wait.");

  if (pos != MPI_UNDEFINED) {
    return pos;
  } else {
    return UInt(-1);
  }
}

/* -------------------------------------------------------------------------- */
void StaticCommunicatorMPI::barrier() {
  MPI_Comm communicator = mpi_data->getMPICommunicator();
  MPI_Barrier(communicator);
}

/* -------------------------------------------------------------------------- */
template <typename T>
void StaticCommunicatorMPI::reduce(T * values, int nb_values,
                                   const SynchronizerOperation & op, int root) {
  MPI_Comm communicator = mpi_data->getMPICommunicator();
  MPI_Datatype type = MPITypeWrapper::getMPIDatatype<T>();

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
      MPI_Reduce(MPI_IN_PLACE, values, nb_values, type,
                 MPITypeWrapper::getMPISynchronizerOperation(op), root,
                 communicator);
  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Allreduce.");
}

/* -------------------------------------------------------------------------- */
template <typename T>
void StaticCommunicatorMPI::allReduce(T * values, int nb_values,
                                      const SynchronizerOperation & op) {
  MPI_Comm communicator = mpi_data->getMPICommunicator();
  MPI_Datatype type = MPITypeWrapper::getMPIDatatype<T>();

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
      MPI_Allreduce(MPI_IN_PLACE, values, nb_values, type,
                    MPITypeWrapper::getMPISynchronizerOperation(op),
                    communicator);
  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Allreduce.");
}

/* -------------------------------------------------------------------------- */
template <typename T>
void StaticCommunicatorMPI::allGather(T * values, int nb_values) {
  MPI_Comm communicator = mpi_data->getMPICommunicator();
  MPI_Datatype type = MPITypeWrapper::getMPIDatatype<T>();

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
      MPI_Allgather(MPI_IN_PLACE, nb_values, type, values, nb_values, type,
                    communicator);

  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Allgather.");
}

/* -------------------------------------------------------------------------- */
template <typename T>
void StaticCommunicatorMPI::allGatherV(T * values, int * nb_values) {
  MPI_Comm communicator = mpi_data->getMPICommunicator();
  int * displs = new int[psize];
  displs[0] = 0;
  for (int i = 1; i < psize; ++i) {
    displs[i] = displs[i - 1] + nb_values[i - 1];
  }

  MPI_Datatype type = MPITypeWrapper::getMPIDatatype<T>();

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
      MPI_Allgatherv(MPI_IN_PLACE, *nb_values, type, values, nb_values, displs,
                     type, communicator);
  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Gather.");

  delete[] displs;
}

/* -------------------------------------------------------------------------- */
template <typename T>
void StaticCommunicatorMPI::gather(T * values, int nb_values, int root) {
  MPI_Comm communicator = mpi_data->getMPICommunicator();
  T *send_buf = NULL, *recv_buf = NULL;
  if (prank == root) {
    send_buf = (T *)MPI_IN_PLACE;
    recv_buf = values;
  } else {
    send_buf = values;
  }

  MPI_Datatype type = MPITypeWrapper::getMPIDatatype<T>();

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
      MPI_Gather(send_buf, nb_values, type, recv_buf, nb_values, type, root,
                 communicator);

  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Gather.");
}

/* -------------------------------------------------------------------------- */
template <typename T>
void StaticCommunicatorMPI::gather(T * values, int nb_values, T * gathered,
                                   int nb_gathered) {
  MPI_Comm communicator = mpi_data->getMPICommunicator();
  T * send_buf = values;
  T * recv_buf = gathered;

  if (nb_gathered == 0)
    nb_gathered = nb_values;

  MPI_Datatype type = MPITypeWrapper::getMPIDatatype<T>();

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
      MPI_Gather(send_buf, nb_values, type, recv_buf, nb_gathered, type,
                 this->prank, communicator);

  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Gather.");
}

/* -------------------------------------------------------------------------- */
template <typename T>
void StaticCommunicatorMPI::gatherV(T * values, int * nb_values, int root) {
  MPI_Comm communicator = mpi_data->getMPICommunicator();
  int * displs = NULL;
  if (prank == root) {
    displs = new int[psize];
    displs[0] = 0;
    for (int i = 1; i < psize; ++i) {
      displs[i] = displs[i - 1] + nb_values[i - 1];
    }
  }

  T *send_buf = NULL, *recv_buf = NULL;
  if (prank == root) {
    send_buf = (T *)MPI_IN_PLACE;
    recv_buf = values;
  } else
    send_buf = values;

  MPI_Datatype type = MPITypeWrapper::getMPIDatatype<T>();

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
      MPI_Gatherv(send_buf, *nb_values, type, recv_buf, nb_values, displs, type,
                  root, communicator);
  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Gather.");

  if (prank == root) {
    delete[] displs;
  }
}

/* -------------------------------------------------------------------------- */
template <typename T>
void StaticCommunicatorMPI::broadcast(T * values, int nb_values, int root) {
  MPI_Comm communicator = mpi_data->getMPICommunicator();
  MPI_Datatype type = MPITypeWrapper::getMPIDatatype<T>();

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
      MPI_Bcast(values, nb_values, type, root, communicator);
  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Gather.");
}

/* -------------------------------------------------------------------------- */
int StaticCommunicatorMPI::getMaxTag() { return this->mpi_data->getMaxTag(); }

/* -------------------------------------------------------------------------- */
int StaticCommunicatorMPI::getMinTag() { return 0; }

/* -------------------------------------------------------------------------- */

// template<typename T>
// MPI_Datatype StaticCommunicatorMPI::getMPIDatatype() {
//   return MPI_DATATYPE_NULL;
// }

template <> MPI_Datatype MPITypeWrapper::getMPIDatatype<char>() {
  return MPI_CHAR;
}

template <> MPI_Datatype MPITypeWrapper::getMPIDatatype<float>() {
  return MPI_FLOAT;
}

template <> MPI_Datatype MPITypeWrapper::getMPIDatatype<double>() {
  return MPI_DOUBLE;
}

template <> MPI_Datatype MPITypeWrapper::getMPIDatatype<long double>() {
  return MPI_LONG_DOUBLE;
}

template <> MPI_Datatype MPITypeWrapper::getMPIDatatype<signed int>() {
  return MPI_INT;
}

template <> MPI_Datatype MPITypeWrapper::getMPIDatatype<NodeType>() {
  return MPITypeWrapper::getMPIDatatype<Int>();
}

template <> MPI_Datatype MPITypeWrapper::getMPIDatatype<unsigned int>() {
  return MPI_UNSIGNED;
}

template <> MPI_Datatype MPITypeWrapper::getMPIDatatype<signed long int>() {
  return MPI_LONG;
}

template <> MPI_Datatype MPITypeWrapper::getMPIDatatype<unsigned long int>() {
  return MPI_UNSIGNED_LONG;
}

template <>
MPI_Datatype MPITypeWrapper::getMPIDatatype<signed long long int>() {
  return MPI_LONG_LONG;
}

template <>
MPI_Datatype MPITypeWrapper::getMPIDatatype<unsigned long long int>() {
  return MPI_UNSIGNED_LONG_LONG;
}

template <>
MPI_Datatype MPITypeWrapper::getMPIDatatype<SCMinMaxLoc<double, int>>() {
  return MPI_DOUBLE_INT;
}

template <>
MPI_Datatype MPITypeWrapper::getMPIDatatype<SCMinMaxLoc<float, int>>() {
  return MPI_FLOAT_INT;
}

template <>
MPI_Datatype MPITypeWrapper::getMPIDatatype<bool>() {
  return MPI_CXX_BOOL;
}

/* -------------------------------------------------------------------------- */
/* Template instantiation                                                     */
/* -------------------------------------------------------------------------- */

#define AKANTU_MPI_COMM_INSTANTIATE(T)                                         \
  template void StaticCommunicatorMPI::send<T>(T * buffer, Int size,           \
                                               Int receiver, Int tag);         \
  template void StaticCommunicatorMPI::receive<T>(T * buffer, Int size,        \
                                                  Int sender, Int tag);        \
  template CommunicationRequest StaticCommunicatorMPI::asyncSend<T>(           \
      T * buffer, Int size, Int receiver, Int tag);                            \
  template CommunicationRequest StaticCommunicatorMPI::asyncReceive<T>(        \
      T * buffer, Int size, Int sender, Int tag);                              \
  template void StaticCommunicatorMPI::probe<T>(Int sender, Int tag,           \
                                                CommunicationStatus & status); \
  template void StaticCommunicatorMPI::allGather<T>(T * values,                \
                                                    int nb_values);            \
  template void StaticCommunicatorMPI::allGatherV<T>(T * values,               \
                                                     int * nb_values);         \
  template void StaticCommunicatorMPI::gather<T>(T * values, int nb_values,    \
                                                 int root);                    \
  template void StaticCommunicatorMPI::gather<T>(                              \
      T * values, int nb_values, T * gathered, int nb_gathered);               \
  template void StaticCommunicatorMPI::gatherV<T>(T * values, int * nb_values, \
                                                  int root);                   \
  template void StaticCommunicatorMPI::broadcast<T>(T * values, int nb_values, \
                                                    int root);                 \
  template void StaticCommunicatorMPI::allReduce<T>(                           \
      T * values, int nb_values, const SynchronizerOperation & op)

AKANTU_MPI_COMM_INSTANTIATE(bool);
AKANTU_MPI_COMM_INSTANTIATE(Real);
AKANTU_MPI_COMM_INSTANTIATE(UInt);
AKANTU_MPI_COMM_INSTANTIATE(Int);
AKANTU_MPI_COMM_INSTANTIATE(char);
AKANTU_MPI_COMM_INSTANTIATE(NodeType);

template void StaticCommunicatorMPI::send<SCMinMaxLoc<Real, int>>(
    SCMinMaxLoc<Real, int> * buffer, Int size, Int receiver, Int tag);
template void StaticCommunicatorMPI::receive<SCMinMaxLoc<Real, int>>(
    SCMinMaxLoc<Real, int> * buffer, Int size, Int sender, Int tag);
template CommunicationRequest
StaticCommunicatorMPI::asyncSend<SCMinMaxLoc<Real, int>>(
    SCMinMaxLoc<Real, int> * buffer, Int size, Int receiver, Int tag);
template CommunicationRequest
StaticCommunicatorMPI::asyncReceive<SCMinMaxLoc<Real, int>>(
    SCMinMaxLoc<Real, int> * buffer, Int size, Int sender, Int tag);
template void StaticCommunicatorMPI::probe<SCMinMaxLoc<Real, int>>(
    Int sender, Int tag, CommunicationStatus & status);
template void StaticCommunicatorMPI::allGather<SCMinMaxLoc<Real, int>>(
    SCMinMaxLoc<Real, int> * values, int nb_values);
template void StaticCommunicatorMPI::allGatherV<SCMinMaxLoc<Real, int>>(
    SCMinMaxLoc<Real, int> * values, int * nb_values);
template void StaticCommunicatorMPI::gather<SCMinMaxLoc<Real, int>>(
    SCMinMaxLoc<Real, int> * values, int nb_values, int root);
template void StaticCommunicatorMPI::gatherV<SCMinMaxLoc<Real, int>>(
    SCMinMaxLoc<Real, int> * values, int * nb_values, int root);
template void StaticCommunicatorMPI::broadcast<SCMinMaxLoc<Real, int>>(
    SCMinMaxLoc<Real, int> * values, int nb_values, int root);
template void StaticCommunicatorMPI::allReduce<SCMinMaxLoc<Real, int>>(
    SCMinMaxLoc<Real, int> * values, int nb_values,
    const SynchronizerOperation & op);

#if AKANTU_INTEGER_SIZE > 4
AKANTU_MPI_COMM_INSTANTIATE(int);
#endif

__END_AKANTU__
