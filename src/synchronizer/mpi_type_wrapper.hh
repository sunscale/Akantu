/**
 * @file   mpi_type_wrapper.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Oct 29 2013
 * @date last modification: Tue Oct 29 2013
 *
 * @brief  Wrapper on MPI types to have a better separation between libraries
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include <mpi.h>
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "static_communicator_mpi.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MPI_TYPE_WRAPPER_HH__
#define __AKANTU_MPI_TYPE_WRAPPER_HH__

__BEGIN_AKANTU__

class MPITypeWrapper {
public:
  MPITypeWrapper(StaticCommunicatorMPI & static_comm) : static_comm(static_comm) {
  }

  template<typename T>
  static inline MPI_Datatype getMPIDatatype();

  inline void setMPICommunicator(MPI_Comm comm) {
    communicator = comm;
    Int prank, psize;
    MPI_Comm_rank(communicator, &prank);
    MPI_Comm_size(communicator, &psize);

    static_comm.setRank(prank);
    static_comm.setSize(psize);
  }

  inline MPI_Comm getMPICommunicator() const {
    return communicator;
  }

  static MPI_Op getMPISynchronizerOperation(const SynchronizerOperation & op) {
    return synchronizer_operation_to_mpi_op[op];
  }

private:
  StaticCommunicatorMPI & static_comm;

  MPI_Comm communicator;

  static MPI_Op synchronizer_operation_to_mpi_op[_so_null + 1];
};

__END_AKANTU__

#endif /* __AKANTU_MPI_TYPE_WRAPPER_HH__ */
