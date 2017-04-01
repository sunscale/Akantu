/**
 * @file   static_communicator_dummy.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Jan 13 2016
 *
 * @brief  Dummy communicator to make everything work im sequential
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

#ifndef __AKANTU_STATIC_COMMUNICATOR_DUMMY_HH__
#define __AKANTU_STATIC_COMMUNICATOR_DUMMY_HH__

/* -------------------------------------------------------------------------- */

#include "aka_common.hh"
#include "real_static_communicator.hh"

/* -------------------------------------------------------------------------- */
#include <cstring>
#include <vector>
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class StaticCommunicatorDummy : public RealStaticCommunicator {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  StaticCommunicatorDummy(int & argc, char **& argv)
      : RealStaticCommunicator(argc, argv) {
    this->prank = 0;
    this->psize = 1;
  }

  virtual ~StaticCommunicatorDummy() {}

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  template <typename T> void send(T *, Int, Int, Int) {}
  template <typename T> void receive(T *, Int, Int, Int) {}

  template <typename T> CommunicationRequest asyncSend(T *, Int, Int, Int) {
    return std::shared_ptr<InternalCommunicationRequest>(
        new InternalCommunicationRequest(0, 0));
  }

  template <typename T> CommunicationRequest asyncReceive(T *, Int, Int, Int) {
    return std::shared_ptr<InternalCommunicationRequest>(
        new InternalCommunicationRequest(0, 0));
  }

  template <typename T> inline void probe(Int, Int, CommunicationStatus &) {}

  bool testRequest(CommunicationRequest &) { return true; }

  void wait(CommunicationRequest &) {}
  void waitAll(std::vector<CommunicationRequest> &) {}
  UInt waitAny(std::vector<CommunicationRequest> &) { return UInt(-1); }

  void barrier() {}

  template <typename T>
  void reduce(T *, int, const SynchronizerOperation &, int) {}

  template <typename T>
  void allReduce(T *, int, const SynchronizerOperation &) {}

  template <typename T> inline void allGather(T *, int) {}
  template <typename T> inline void allGatherV(T *, int *) {}

  template <typename T> inline void gather(T *, int, int = 0) {}
  template <typename T>
  inline void gather(T * values, int nb_values, T * gathered, int) {
    std::memcpy(gathered, values, nb_values);
  }

  template <typename T> inline void gatherV(T *, int *, int = 0) {}
  template <typename T> inline void broadcast(T *, int, int = 0) {}

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  int getMaxTag() { return std::numeric_limits<int>::max(); }
  int getMinTag() { return 0; }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
};

__END_AKANTU__

#endif /* __AKANTU_STATIC_COMMUNICATOR_DUMMY_HH__ */
