/**
 * @file   static_communicator.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Jan 19 2016
 *
 * @brief  Class handling the parallel communications
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

#ifndef __AKANTU_STATIC_COMMUNICATOR_HH__
#define __AKANTU_STATIC_COMMUNICATOR_HH__

/* -------------------------------------------------------------------------- */
#include "aka_array.hh"
#include "aka_common.hh"
#include "aka_event_handler_manager.hh"
#include "communication_buffer.hh"
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
#define AKANTU_COMMUNICATOR_LIST_0 BOOST_PP_SEQ_NIL

#include "static_communicator_dummy.hh"
#define AKANTU_COMMUNICATOR_LIST_1                                             \
  BOOST_PP_SEQ_PUSH_BACK(                                                      \
      AKANTU_COMMUNICATOR_LIST_0,                                              \
      (_communicator_dummy, (StaticCommunicatorDummy, BOOST_PP_NIL)))

#if defined(AKANTU_USE_MPI)
#include "static_communicator_mpi.hh"
#define AKANTU_COMMUNICATOR_LIST_ALL                                           \
  BOOST_PP_SEQ_PUSH_BACK(                                                      \
      AKANTU_COMMUNICATOR_LIST_1,                                              \
      (_communicator_mpi, (StaticCommunicatorMPI, BOOST_PP_NIL)))
#else
#define AKANTU_COMMUNICATOR_LIST_ALL AKANTU_COMMUNICATOR_LIST_1
#endif // AKANTU_COMMUNICATOR_LIST

#include "real_static_communicator.hh"

/* -------------------------------------------------------------------------- */

namespace akantu {

class RealStaticCommunicator;

struct FinalizeCommunicatorEvent {
  explicit FinalizeCommunicatorEvent(const StaticCommunicator & comm)
      : communicator(comm) {}
  const StaticCommunicator & communicator;
};

class CommunicatorEventHandler {
public:
  virtual ~CommunicatorEventHandler() {}
  virtual void onCommunicatorFinalize() {}

private:
  inline void sendEvent(const FinalizeCommunicatorEvent &) {
    onCommunicatorFinalize();
  }

  template <class EventHandler> friend class EventHandlerManager;
};

class StaticCommunicator
    : public EventHandlerManager<CommunicatorEventHandler> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
protected:
  StaticCommunicator(int & argc, char **& argv,
                     CommunicatorType type = _communicator_mpi);

public:
  virtual ~StaticCommunicator();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Point to Point                                                           */
  /* ------------------------------------------------------------------------ */
  template <typename T>
  inline void receive(Array<T> & values, Int sender, Int tag) const {
    return this->receive(values.storage(),
                         values.getSize() * values.getNbComponent(), sender,
                         tag);
  }
  template <typename T>
  inline void receive(Vector<T> & values, Int sender, Int tag) const {
    return this->receive(values.storage(), values.size(), sender, tag);
  }
  template <typename T>
  inline void receive(Matrix<T> & values, Int sender, Int tag) const {
    return this->receive(values.storage(), values.size(), sender, tag);
  }
  template <bool is_static>
  inline void receive(CommunicationBufferTemplated<is_static> & values,
                      Int sender, Int tag) const {
    return this->receive(values.storage(), values.getSize(), sender, tag);
  }
  template <typename T>
  inline void receive(T & values, Int sender, Int tag) const {
    return this->receive(&values, 1, sender, tag);
  }
  /* ------------------------------------------------------------------------ */
  template <typename T>
  inline void send(Array<T> & values, Int receiver, Int tag) const {
    return this->send(values.storage(),
                      values.getSize() * values.getNbComponent(), receiver,
                      tag);
  }
  template <typename T>
  inline void send(Vector<T> & values, Int receiver, Int tag) const {
    return this->send(values.storage(), values.size(), receiver, tag);
  }
  template <typename T>
  inline void send(Matrix<T> & values, Int receiver, Int tag) const {
    return this->send(values.storage(), values.size(), receiver, tag);
  }
  template <bool is_static>
  inline void send(CommunicationBufferTemplated<is_static> & values,
                   Int receiver, Int tag) const {
    return this->send(values.storage(), values.getSize(), receiver, tag);
  }
  template <typename T>
  inline void send(T & values, Int receiver, Int tag) const {
    return this->send(&values, 1, receiver, tag);
  }

  /* ------------------------------------------------------------------------ */
  template <typename T>
  inline CommunicationRequest asyncSend(Array<T> & values, Int receiver,
                                        Int tag) const {
    return this->asyncSend(values.storage(),
                           values.getSize() * values.getNbComponent(), receiver,
                           tag);
  }
  template <typename T>
  inline CommunicationRequest asyncSend(std::vector<T> & values, Int receiver,
                                        Int tag) const {
    return this->asyncSend(values.data(), values.size(), receiver, tag);
  }

  template <typename T>
  inline CommunicationRequest asyncSend(Vector<T> & values, Int receiver,
                                        Int tag) const {
    return this->asyncSend(values.storage(), values.size(), receiver, tag);
  }
  template <typename T>
  inline CommunicationRequest asyncSend(Matrix<T> & values, Int receiver,
                                        Int tag) const {
    return this->asyncSend(values.storage(), values.size(), receiver, tag);
  }
  template <bool is_static>
  inline CommunicationRequest
  asyncSend(CommunicationBufferTemplated<is_static> & values, Int receiver,
            Int tag) const {
    return this->asyncSend(values.storage(), values.getSize(), receiver, tag);
  }
  template <typename T>
  inline CommunicationRequest asyncSend(T & values, Int receiver,
                                        Int tag) const {
    return this->asyncSend(&values, 1, receiver, tag);
  }

  /* ------------------------------------------------------------------------ */
  template <typename T>
  inline CommunicationRequest asyncReceive(Array<T> & values, Int sender,
                                           Int tag) const {
    return this->asyncReceive(values.storage(),
                              values.getSize() * values.getNbComponent(),
                              sender, tag);
  }
  template <typename T, UInt ndim, class RetType>
  inline CommunicationRequest
  asyncReceive(TensorStorage<T, ndim, RetType> & values, Int sender,
               Int tag) const {
    return this->asyncReceive(values.storage(), values.size(), sender, tag);
  }
  template <bool is_static>
  inline CommunicationRequest
  asyncReceive(CommunicationBufferTemplated<is_static> & values, Int sender,
               Int tag) const {
    return this->asyncReceive(values.storage(), values.getSize(), sender, tag);
  }

  /* ------------------------------------------------------------------------ */
  template <typename T>
  inline void probe(Int sender, Int tag, CommunicationStatus & status) const;

  /* ------------------------------------------------------------------------ */
  /* Collectives                                                              */
  /* ------------------------------------------------------------------------ */
  template <typename T>
  inline void allReduce(Array<T> & values,
                        const SynchronizerOperation & op) const {
    this->allReduce(values.storage(),
                    values.getSize() * values.getNbComponent(), op);
  }
  template <typename T>
  inline void allReduce(Vector<T> & values,
                        const SynchronizerOperation & op) const {
    this->allReduce(values.storage(), values.size(), op);
  }
  template <typename T>
  inline void allReduce(Matrix<T> & values,
                        const SynchronizerOperation & op) const {
    this->allReduce(values.storage(), values.size(), op);
  }
  template <typename T>
  inline void allReduce(T & values, const SynchronizerOperation & op) const {
    this->allReduce(&values, 1, op);
  }

  /* ------------------------------------------------------------------------ */
  template <typename T> inline void allGather(Array<T> & values) const {
    AKANTU_DEBUG_ASSERT(UInt(this->real_static_communicator->getNbProc()) ==
                            values.getSize(),
                        "The array size is not correct");
    this->allGather(values.storage(), values.getNbComponent());
  }

  template <typename T> inline void allGather(Vector<T> & values) const {
    AKANTU_DEBUG_ASSERT(UInt(this->real_static_communicator->getNbProc()) ==
                            values.size(),
                        "The vector size is not correct");
    this->allGather(values.storage(), 1);
  }

  template <typename T> inline void allGather(Tensor3<T> & values) const {
    AKANTU_DEBUG_ASSERT(UInt(this->real_static_communicator->getNbProc()) ==
                            values.size(2),
                        "The tensor size is not correct");
    this->allGather(values.storage(), values.size(0) * values.size(1));
  }

  /* ------------------------------------------------------------------------ */
  template <typename T>
  inline void allGatherV(Array<T> & values, Array<Int> sizes) const {
    this->allGatherV(values.storage(), sizes.storage());
  }

  /* ------------------------------------------------------------------------ */
  template <typename T>
  inline void reduce(Array<T> & values, const SynchronizerOperation & op,
                     int root = 0) const {
    this->reduce(values.storage(), values.getSize() * values.getNbComponent(),
                 op, root);
  }

  /* ------------------------------------------------------------------------ */
  template <typename T>
  inline void gather(Vector<T> & values, int root = 0) const {
    this->gather(values.storage(), values.getNbComponent(), root);
  }
  template <typename T> inline void gather(T values, int root = 0) const {
    this->gather(&values, 1, root);
  }
  /* ------------------------------------------------------------------------ */
  template <typename T>
  inline void gather(Vector<T> & values, Array<T> & gathered) const {
    AKANTU_DEBUG_ASSERT(values.getSize() == gathered.getNbComponent(),
                        "The array size is not correct");
    gathered.resize(this->real_static_communicator->getNbProc());
    this->gather(values.data(), values.getSize(), gathered.storage(),
                 gathered.getNbComponent());
  }

  template <typename T>
  inline void gather(T values, Array<T> & gathered) const {
    this->gather(&values, 1, gathered.storage(), 1);
  }

  /* ------------------------------------------------------------------------ */
  template <typename T>
  inline void gatherV(Array<T> & values, Array<Int> sizes, int root = 0) const {
    this->gatherV(values.storage(), sizes.storage(), root);
  }

  /* ------------------------------------------------------------------------ */
  template <typename T>
  inline void broadcast(Array<T> & values, int root = 0) const {
    this->broadcast(values.storage(),
                    values.getSize() * values.getNbComponent(), root);
  }
  template <bool is_static>
  inline void broadcast(CommunicationBufferTemplated<is_static> & values,
                        int root = 0) const {
    this->broadcast(values.storage(), values.getSize(), root);
  }
  template <typename T> inline void broadcast(T & values, int root = 0) const {
    this->broadcast(&values, 1, root);
  }

  /* ------------------------------------------------------------------------ */
  inline void barrier() const;

  /* ------------------------------------------------------------------------ */
  /* Request handling                                                         */
  /* ------------------------------------------------------------------------ */
  inline bool testRequest(CommunicationRequest & request) const;

  inline void wait(CommunicationRequest & request) const;
  inline void waitAll(std::vector<CommunicationRequest> & requests) const;
  inline UInt waitAny(std::vector<CommunicationRequest> & requests) const;

  inline void freeCommunicationRequest(CommunicationRequest & request) const;
  inline void
  freeCommunicationRequest(std::vector<CommunicationRequest> & requests) const;

protected:
  template <typename T>
  inline void send(T * buffer, Int size, Int receiver, Int tag) const;
  template <typename T>
  inline void receive(T * buffer, Int size, Int sender, Int tag) const;

  template <typename T>
  inline CommunicationRequest asyncSend(T * buffer, Int size, Int receiver,
                                        Int tag) const;
  template <typename T>
  inline CommunicationRequest asyncReceive(T * buffer, Int size, Int sender,
                                           Int tag) const;

  template <typename T>
  inline void allReduce(T * values, int nb_values,
                        const SynchronizerOperation & op) const;

  template <typename T> inline void allGather(T * values, int nb_values) const;
  template <typename T>
  inline void allGatherV(T * values, int * nb_values) const;

  template <typename T>
  inline void reduce(T * values, int nb_values,
                     const SynchronizerOperation & op, int root = 0) const;
  template <typename T>
  inline void gather(T * values, int nb_values, int root = 0) const;

  template <typename T>
  inline void gather(T * values, int nb_values, T * gathered,
                     int nb_gathered = 0) const;

  template <typename T>
  inline void gatherV(T * values, int * nb_values, int root = 0) const;

  template <typename T>
  inline void broadcast(T * values, int nb_values, int root = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  virtual Int getNbProc() const { return real_static_communicator->psize; };
  virtual Int whoAmI() const { return real_static_communicator->prank; };

  AKANTU_GET_MACRO(RealStaticCommunicator, *real_static_communicator,
                   const RealStaticCommunicator &);
  AKANTU_GET_MACRO_NOT_CONST(RealStaticCommunicator, *real_static_communicator,
                             RealStaticCommunicator &);

  template <class Comm> Comm & getRealStaticCommunicator() {
    return dynamic_cast<Comm &>(*real_static_communicator);
  }

  static StaticCommunicator &
  getStaticCommunicator(CommunicatorType type = _communicator_mpi);

  static StaticCommunicator &
  getStaticCommunicator(int & argc, char **& argv,
                        CommunicatorType type = _communicator_mpi);

  static StaticCommunicator * getStaticCommunicatorDummy();

  static bool isInstantiated() { return is_instantiated; };

  int getMaxTag() const;
  int getMinTag() const;
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  static bool is_instantiated;

  static StaticCommunicator * static_communicator;

  RealStaticCommunicator * real_static_communicator;

  CommunicatorType real_type;
};

/* -------------------------------------------------------------------------- */
/* Inline Functions ArrayBase                                                */
/* -------------------------------------------------------------------------- */
inline std::ostream & operator<<(std::ostream & stream,
                                 const CommunicationRequest & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#include "static_communicator_inline_impl.hh"

#endif /* __AKANTU_STATIC_COMMUNICATOR_HH__ */
