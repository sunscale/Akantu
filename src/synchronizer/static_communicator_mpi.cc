/**
 * @file   static_communicator_mpi.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sun Sep 26 2010
 * @date last modification: Mon Jul 21 2014
 *
 * @brief  StaticCommunicatorMPI implementation
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include <iostream>

/* -------------------------------------------------------------------------- */
#include "static_communicator_mpi.hh"
#include "mpi_type_wrapper.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

MPI_Op MPITypeWrapper::synchronizer_operation_to_mpi_op[_so_null + 1] = {
  MPI_SUM,
  MPI_MIN,
  MPI_MAX,
  MPI_OP_NULL
};


class CommunicationRequestMPI : public CommunicationRequest {
public:
  CommunicationRequestMPI(UInt source, UInt dest);
  ~CommunicationRequestMPI();
  MPI_Request * getMPIRequest() { return request; };
private:
  MPI_Request * request;
};


/* -------------------------------------------------------------------------- */
/* Implementation                                                             */
/* -------------------------------------------------------------------------- */
CommunicationRequestMPI::CommunicationRequestMPI(UInt source, UInt dest) :
  CommunicationRequest(source, dest) {
  request = new MPI_Request;
}

/* -------------------------------------------------------------------------- */
CommunicationRequestMPI::~CommunicationRequestMPI() {
  delete request;
}

/* -------------------------------------------------------------------------- */
StaticCommunicatorMPI::StaticCommunicatorMPI(int & argc, char ** & argv) :
  RealStaticCommunicator(argc, argv) {

  if(argc != 0) {
    MPI_Init(&argc, &argv);
  }

  mpi_data = new MPITypeWrapper(*this);
  mpi_data->setMPICommunicator(MPI_COMM_WORLD);
}

/* -------------------------------------------------------------------------- */
StaticCommunicatorMPI::~StaticCommunicatorMPI() {
  MPI_Finalize();
}

/* -------------------------------------------------------------------------- */
template<typename T>
void StaticCommunicatorMPI::send(T * buffer, Int size,
				 Int receiver, Int tag) {
  MPI_Comm communicator = mpi_data->getMPICommunicator();
  MPI_Datatype type = MPITypeWrapper::getMPIDatatype<T>();
#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
    MPI_Send(buffer, size, type, receiver, tag, communicator);

  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Send.");
}

/* -------------------------------------------------------------------------- */
template<typename T>
void StaticCommunicatorMPI::receive(T * buffer, Int size,
				    Int sender, Int tag) {
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
template<typename T>
CommunicationRequest * StaticCommunicatorMPI::asyncSend(T * buffer, Int size,
							Int receiver, Int tag) {
  MPI_Comm communicator = mpi_data->getMPICommunicator();
  CommunicationRequestMPI * request = new CommunicationRequestMPI(prank, receiver);
  MPI_Datatype type = MPITypeWrapper::getMPIDatatype<T>();

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
    MPI_Isend(buffer, size, type, receiver, tag, communicator, request->getMPIRequest());

  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Isend.");
  return request;
}

/* -------------------------------------------------------------------------- */
template<typename T>
CommunicationRequest * StaticCommunicatorMPI::asyncReceive(T * buffer, Int size,
							   Int sender, Int tag) {
  MPI_Comm communicator = mpi_data->getMPICommunicator();
  CommunicationRequestMPI * request = new CommunicationRequestMPI(sender, prank);
  MPI_Datatype type = MPITypeWrapper::getMPIDatatype<T>();

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
    MPI_Irecv(buffer, size, type, sender, tag, communicator, request->getMPIRequest());

  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Irecv.");
  return request;
}

/* -------------------------------------------------------------------------- */
template<typename T>
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
bool StaticCommunicatorMPI::testRequest(CommunicationRequest * request) {
  MPI_Status status;
  int flag;
  CommunicationRequestMPI * req_mpi = static_cast<CommunicationRequestMPI *>(request);
  MPI_Request * req = req_mpi->getMPIRequest();

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif

    MPI_Test(req, &flag, &status);

  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Test.");
  return (flag != 0);
}

/* -------------------------------------------------------------------------- */
void StaticCommunicatorMPI::wait(CommunicationRequest * request) {
  MPI_Status status;
  CommunicationRequestMPI * req_mpi = static_cast<CommunicationRequestMPI *>(request);
  MPI_Request * req = req_mpi->getMPIRequest();

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
    MPI_Wait(req, &status);

  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Wait.");
}

/* -------------------------------------------------------------------------- */
void StaticCommunicatorMPI::waitAll(std::vector<CommunicationRequest *> & requests) {
  MPI_Status status;
  std::vector<CommunicationRequest *>::iterator it;
  for(it = requests.begin(); it != requests.end(); ++it) {
    MPI_Request * req = static_cast<CommunicationRequestMPI *>(*it)->getMPIRequest();

#if !defined(AKANTU_NDEBUG)
    int ret =
#endif
      MPI_Wait(req, &status);

    AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Wait.");
  }
}

/* -------------------------------------------------------------------------- */
void StaticCommunicatorMPI::barrier() {
  MPI_Comm communicator = mpi_data->getMPICommunicator();
  MPI_Barrier(communicator);
}

/* -------------------------------------------------------------------------- */
template<typename T>
void StaticCommunicatorMPI::allReduce(T * values, Int nb_values,
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
template<typename T>
void StaticCommunicatorMPI::allGather(T * values, Int nb_values) {
  MPI_Comm communicator = mpi_data->getMPICommunicator();
  MPI_Datatype type = MPITypeWrapper::getMPIDatatype<T>();

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
    MPI_Allgather(MPI_IN_PLACE, nb_values, type, values, nb_values, type, communicator);

  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Allgather.");
}

/* -------------------------------------------------------------------------- */
template<typename T>
void StaticCommunicatorMPI::allGatherV(T * values, Int * nb_values) {
  MPI_Comm communicator = mpi_data->getMPICommunicator();
  Int * displs = new Int[psize];
  displs[0] = 0;
  for (Int i = 1; i < psize; ++i) {
    displs[i] = displs[i-1] + nb_values[i-1];
  }

  MPI_Datatype type = MPITypeWrapper::getMPIDatatype<T>();

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
    MPI_Allgatherv(MPI_IN_PLACE, *nb_values, type, values, nb_values, displs, type, communicator);
  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Gather.");

  delete [] displs;
}

/* -------------------------------------------------------------------------- */
template<typename T>
void StaticCommunicatorMPI::gather(T * values, Int nb_values, Int root) {
  MPI_Comm communicator = mpi_data->getMPICommunicator();
  T * send_buf = NULL, * recv_buf = NULL;
  if(prank == root) {
    send_buf = (T *) MPI_IN_PLACE;
    recv_buf = values;
  } else {
    send_buf = values;
  }

  MPI_Datatype type = MPITypeWrapper::getMPIDatatype<T>();

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
    MPI_Gather(send_buf, nb_values, type, recv_buf, nb_values, type, root, communicator);

  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Gather.");
}

/* -------------------------------------------------------------------------- */
template<typename T>
void StaticCommunicatorMPI::gatherV(T * values, Int * nb_values, Int root) {
  MPI_Comm communicator = mpi_data->getMPICommunicator();
  Int * displs = NULL;
  if(prank == root) {
    displs = new Int[psize];
    displs[0] = 0;
    for (Int i = 1; i < psize; ++i) {
      displs[i] = displs[i-1] + nb_values[i-1];
    }
  }

  T * send_buf = NULL, * recv_buf = NULL;
  if(prank == root) {
    send_buf = (T *) MPI_IN_PLACE;
    recv_buf = values;
  } else send_buf = values;

  MPI_Datatype type = MPITypeWrapper::getMPIDatatype<T>();

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
    MPI_Gatherv(send_buf, *nb_values, type, recv_buf, nb_values, displs, type, root, communicator);
  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Gather.");

  if(prank == root) {
    delete [] displs;
  }
}

/* -------------------------------------------------------------------------- */
template<typename T>
void StaticCommunicatorMPI::broadcast(T * values, Int nb_values, Int root) {
  MPI_Comm communicator = mpi_data->getMPICommunicator();
  MPI_Datatype type = MPITypeWrapper::getMPIDatatype<T>();

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
    MPI_Bcast(values, nb_values, type, root, communicator);
  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Gather.");
}

/* -------------------------------------------------------------------------- */
// template<typename T>
// MPI_Datatype StaticCommunicatorMPI::getMPIDatatype() {
//   return MPI_DATATYPE_NULL;
// }

template<>
MPI_Datatype MPITypeWrapper::getMPIDatatype<char>() {
  return MPI_CHAR;
}

template<>
MPI_Datatype MPITypeWrapper::getMPIDatatype<Real>() {
  return MPI_DOUBLE;
}

template<>
MPI_Datatype MPITypeWrapper::getMPIDatatype<UInt>() {
  return MPI_UNSIGNED;
}

template<>
MPI_Datatype MPITypeWrapper::getMPIDatatype<Int>() {
  return MPI_INT;
}

/* -------------------------------------------------------------------------- */
/* Template instantiation                                                     */
/* -------------------------------------------------------------------------- */

#define AKANTU_MPI_COMM_INSTANTIATE(T)					\
  template void StaticCommunicatorMPI::send<T>   (T * buffer, Int size, Int receiver, Int tag); \
  template void StaticCommunicatorMPI::receive<T>(T * buffer, Int size, Int sender,   Int tag); \
  template CommunicationRequest * StaticCommunicatorMPI::asyncSend<T>   (T * buffer, Int size, Int receiver, Int tag); \
  template CommunicationRequest * StaticCommunicatorMPI::asyncReceive<T>(T * buffer, Int size, Int sender,   Int tag); \
  template void StaticCommunicatorMPI::probe<T>(Int sender, Int tag, CommunicationStatus & status); \
  template void StaticCommunicatorMPI::allGather<T> (T * values, Int nb_values); \
  template void StaticCommunicatorMPI::allGatherV<T>(T * values, Int * nb_values); \
  template void StaticCommunicatorMPI::gather<T> (T * values, Int nb_values, Int root);	\
  template void StaticCommunicatorMPI::gatherV<T>(T * values, Int * nb_values, Int root); \
  template void StaticCommunicatorMPI::broadcast<T>(T * values, Int nb_values, Int root); \
  template void StaticCommunicatorMPI::allReduce<T>(T * values, Int nb_values, const SynchronizerOperation & op);


AKANTU_MPI_COMM_INSTANTIATE(Real);
AKANTU_MPI_COMM_INSTANTIATE(UInt);
AKANTU_MPI_COMM_INSTANTIATE(Int);
AKANTU_MPI_COMM_INSTANTIATE(char);


__END_AKANTU__
