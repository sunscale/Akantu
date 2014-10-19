/**
 * @file   static_communicator_mpi.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sun Sep 05 2010
 * @date last modification: Tue Oct 29 2013
 *
 * @brief  class handling parallel communication trough MPI
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

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_STATIC_COMMUNICATOR_MPI_HH__
#define __AKANTU_STATIC_COMMUNICATOR_MPI_HH__

/* -------------------------------------------------------------------------- */
#include "real_static_communicator.hh"

/* -------------------------------------------------------------------------- */
#include <vector>


__BEGIN_AKANTU__

class MPITypeWrapper;

/* -------------------------------------------------------------------------- */
class StaticCommunicatorMPI : public RealStaticCommunicator {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  StaticCommunicatorMPI(int & argc, char ** & argv);

  virtual ~StaticCommunicatorMPI();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  template<typename T> void send   (T * buffer, Int size, Int receiver, Int tag);
  template<typename T> void receive(T * buffer, Int size, Int sender,   Int tag);

  template<typename T> CommunicationRequest * asyncSend   (T * buffer,   Int size,
							   Int receiver, Int tag);
  template<typename T> CommunicationRequest * asyncReceive(T * buffer,   Int size,
							   Int sender,   Int tag);

  template<typename T> void probe(Int sender, Int tag,
				  CommunicationStatus & status);

  template<typename T> void allGather (T * values, Int nb_values);
  template<typename T> void allGatherV(T * values, Int * nb_values);

  template<typename T> void gather (T * values, Int nb_values, Int root);
  template<typename T> void gatherV(T * values, Int * nb_values, Int root);

  template<typename T> void broadcast(T * values, Int nb_values, Int root);

  bool testRequest(CommunicationRequest * request);

  void wait   (CommunicationRequest * request);
  void waitAll(std::vector<CommunicationRequest *> & requests);

  void barrier();

  template<typename T> void allReduce(T * values, Int nb_values,
				      const SynchronizerOperation & op);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  const MPITypeWrapper & getMPITypeWrapper() const { return *mpi_data; }
  MPITypeWrapper & getMPITypeWrapper() { return *mpi_data; }

private:
  void setRank(int prank) { this->prank = prank; }
  void setSize(int psize) { this->psize = psize; }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  friend class MPITypeWrapper;

  MPITypeWrapper * mpi_data;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "static_communicator_mpi_inline_impl.hh"

__END_AKANTU__

#endif /* __AKANTU_STATIC_COMMUNICATOR_MPI_HH__ */
