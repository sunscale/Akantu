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
#include <vector>
#include <cstring>
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class StaticCommunicatorDummy : public RealStaticCommunicator {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  StaticCommunicatorDummy(__attribute__((unused)) int & argc,
                          __attribute__((unused)) char **& argv)
      : RealStaticCommunicator(argc, argv) {
    this->prank = 0;
    this->psize = 1;
  }

  virtual ~StaticCommunicatorDummy() {}

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  template <typename T>
  void
  send(__attribute__((unused)) T * buffer, __attribute__((unused)) Int size,
       __attribute__((unused)) Int receiver, __attribute__((unused)) Int tag) {}

  template <typename T>
  void receive(__attribute__((unused)) T * buffer,
               __attribute__((unused)) Int size,
               __attribute__((unused)) Int sender,
               __attribute__((unused)) Int tag) {}

  template <typename T>
  CommunicationRequest asyncSend(__attribute__((unused)) T * buffer,
                                   __attribute__((unused)) Int size,
                                   __attribute__((unused)) Int receiver,
                                   __attribute__((unused)) Int tag) {
    return (new InternalCommunicationRequest(0, 0));
  }

  template <typename T>
  CommunicationRequest asyncReceive(__attribute__((unused)) T * buffer,
                                      __attribute__((unused)) Int size,
                                      __attribute__((unused)) Int sender,
                                      __attribute__((unused)) Int tag) {
    return (new InternalCommunicationRequest(0, 0));
  }

  template <typename T>
  inline void probe(__attribute__((unused)) Int sender,
                    __attribute__((unused)) Int tag,
                    __attribute__((unused)) CommunicationStatus & status) {}

  bool testRequest(__attribute__((unused)) CommunicationRequest request) {
    return true;
  }

  void wait(__attribute__((unused)) CommunicationRequest request) {}
  void waitAll(__attribute__((unused))
               std::vector<CommunicationRequest> & requests) {}
  UInt waitAny(__attribute__((unused))
               std::vector<CommunicationRequest> & requests) { return UInt(-1); }


  void barrier() {}

  template <typename T>
  void reduce(__attribute__((unused)) T * values,
              __attribute__((unused)) int nb_values,
              __attribute__((unused)) const SynchronizerOperation & op,
              __attribute__((unused)) int root) {}

  template <typename T>
  void allReduce(__attribute__((unused)) T * values,
                 __attribute__((unused)) int nb_values,
                 __attribute__((unused)) const SynchronizerOperation & op) {}

  template <typename T>
  inline void allGather(__attribute__((unused)) T * values,
                        __attribute__((unused)) int nb_values) {}

  template <typename T>
  inline void allGatherV(__attribute__((unused)) T * values,
                         __attribute__((unused)) int * nb_values) {}

  template <typename T>
  inline void gather(__attribute__((unused)) T * values,
                     __attribute__((unused)) int nb_values,
                     __attribute__((unused)) int root = 0) {}

  template <typename T>
  inline void gather(T * values,
                     int nb_values,
                     T * gathered,
                     __attribute__((unused)) int nb_gathered) {
    std::memcpy(gathered, values, nb_values);
  }

  template <typename T>
  inline void gatherV(__attribute__((unused)) T * values,
                      __attribute__((unused)) int * nb_values,
                      __attribute__((unused)) int root = 0) {}

  template <typename T>
  inline void broadcast(__attribute__((unused)) T * values,
                        __attribute__((unused)) int nb_values,
                        __attribute__((unused)) int root = 0) {}

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
