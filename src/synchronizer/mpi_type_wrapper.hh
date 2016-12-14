/**
 * @file   mpi_type_wrapper.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Jun 14 2010
 * @date last modification: Wed Oct 07 2015
 *
 * @brief  Wrapper on MPI types to have a better separation between libraries
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
#if defined(__INTEL_COMPILER)
//#pragma warning ( disable : 383 )
#elif defined(__clang__) // test clang to be sure that when we test for gnu it
                         // is only gnu
#elif (defined(__GNUC__) || defined(__GNUG__))
#if __cplusplus > 199711L
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wliteral-suffix"
#endif
#endif

#include <mpi.h>

#if defined(__INTEL_COMPILER)
//#pragma warning ( disable : 383 )
#elif defined(__clang__) // test clang to be sure that when we test for gnu it
                         // is only gnu
#elif (defined(__GNUC__) || defined(__GNUG__))
#if __cplusplus > 199711L
#pragma GCC diagnostic pop
#endif
#endif

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "static_communicator_mpi.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MPI_TYPE_WRAPPER_HH__
#define __AKANTU_MPI_TYPE_WRAPPER_HH__

__BEGIN_AKANTU__

class MPITypeWrapper {
public:
  MPITypeWrapper(StaticCommunicatorMPI & static_comm)
      : static_comm(static_comm), communicator(MPI_COMM_WORLD), max_tag(0) {}

  template <typename T> static inline MPI_Datatype getMPIDatatype();

  inline void setMPICommunicator(MPI_Comm comm) {
    communicator = comm;
    int prank, psize;
    MPI_Comm_rank(communicator, &prank);
    MPI_Comm_size(communicator, &psize);

    static_comm.setRank(prank);
    static_comm.setSize(psize);

    int flag;
    int * value;
    MPI_Comm_get_attr(comm, MPI_TAG_UB, &value, &flag);
    AKANTU_DEBUG_ASSERT(flag, "No attribute MPI_TAG_UB.");
    this->max_tag = *value;
  }

  inline MPI_Comm getMPICommunicator() const { return communicator; }

  static MPI_Op getMPISynchronizerOperation(const SynchronizerOperation & op) {
    return synchronizer_operation_to_mpi_op[op];
  }

  inline int getMaxTag() const { return this->max_tag; }

private:
  StaticCommunicatorMPI & static_comm;

  MPI_Comm communicator;

  int max_tag;

  static MPI_Op synchronizer_operation_to_mpi_op[_so_null + 1];
};

__END_AKANTU__

#endif /* __AKANTU_MPI_TYPE_WRAPPER_HH__ */
