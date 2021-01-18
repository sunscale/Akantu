/**
 * @file   mpi_communicator_data.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Jun 14 2010
 * @date last modification: Mon Feb 05 2018
 *
 * @brief  Wrapper on MPI types to have a better separation between libraries
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "communicator.hh"
/* -------------------------------------------------------------------------- */
#include <unordered_map>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MPI_TYPE_WRAPPER_HH_
#define AKANTU_MPI_TYPE_WRAPPER_HH_

namespace akantu {

class MPICommunicatorData : public CommunicatorInternalData {
public:
  MPICommunicatorData() {
    MPI_Initialized(&is_externaly_initialized);
    if (is_externaly_initialized == 0) {
      MPI_Init(nullptr, nullptr); // valid according to the spec
    }

    MPI_Comm_create_errhandler(MPICommunicatorData::errorHandler,
                               &error_handler);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, error_handler);
    setMPICommunicator(MPI_COMM_WORLD);
  }

  ~MPICommunicatorData() override {
    if (is_externaly_initialized == 0) {
      MPI_Comm_set_errhandler(communicator, save_error_handler);
      MPI_Errhandler_free(&error_handler);
      MPI_Finalize();
    }
  }

  inline void setMPICommunicator(MPI_Comm comm) {
    MPI_Comm_set_errhandler(communicator, save_error_handler);
    communicator = comm;
    MPI_Comm_get_errhandler(comm, &save_error_handler);
    MPI_Comm_set_errhandler(comm, error_handler);
  }

  inline int rank() const {
    int prank;
    MPI_Comm_rank(communicator, &prank);
    return prank;
  }

  inline int size() const {
    int psize;
    MPI_Comm_size(communicator, &psize);
    return psize;
  }

  inline MPI_Comm getMPICommunicator() const { return communicator; }
  static int getMaxTag() {
    int flag;
    int * value;
    // not defined on derived intra-communicator
    MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &value, &flag);
    AKANTU_DEBUG_ASSERT(flag, "No attribute MPI_TAG_UB.");
    return *value;
  }

private:
  MPI_Comm communicator{MPI_COMM_WORLD};
  MPI_Errhandler save_error_handler{MPI_ERRORS_ARE_FATAL};
  static int is_externaly_initialized;
  /* ------------------------------------------------------------------------ */
  MPI_Errhandler error_handler;

  static void
  errorHandler(MPI_Comm * /*comm*/,
               int * error_code, // NOLINT(readability-non-const-parameter)
               ...) {
    char error_string[MPI_MAX_ERROR_STRING];
    int str_len;
    MPI_Error_string(*error_code, error_string, &str_len);
    AKANTU_CUSTOM_EXCEPTION_INFO(debug::CommunicationException(),
                                 "MPI failed with the error code "
                                     << *error_code << ": \"" << error_string
                                     << "\"");
  }
};

} // namespace akantu

#endif /* AKANTU_MPI_TYPE_WRAPPER_HH_ */
