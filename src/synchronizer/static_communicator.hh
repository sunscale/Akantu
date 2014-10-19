/**
 * @file   static_communicator.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Sep 01 2010
 * @date last modification: Mon Jul 21 2014
 *
 * @brief  Class handling the parallel communications
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

#ifndef __AKANTU_STATIC_COMMUNICATOR_HH__
#define __AKANTU_STATIC_COMMUNICATOR_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_event_handler_manager.hh"

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
#define AKANTU_COMMUNICATOR_LIST_0 BOOST_PP_SEQ_NIL

#include "static_communicator_dummy.hh"
#define AKANTU_COMMUNICATOR_LIST_1					\
  BOOST_PP_SEQ_PUSH_BACK(AKANTU_COMMUNICATOR_LIST_0,			\
			 (_communicator_dummy,	(StaticCommunicatorDummy, BOOST_PP_NIL)))

#if defined(AKANTU_USE_MPI)
#  include "static_communicator_mpi.hh"
#  define AKANTU_COMMUNICATOR_LIST_ALL					\
  BOOST_PP_SEQ_PUSH_BACK(AKANTU_COMMUNICATOR_LIST_1,			\
			 (_communicator_mpi, (StaticCommunicatorMPI, BOOST_PP_NIL)))
#else
#  define AKANTU_COMMUNICATOR_LIST_ALL	AKANTU_COMMUNICATOR_LIST_1
#endif // AKANTU_COMMUNICATOR_LIST

#include "real_static_communicator.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class RealStaticCommunicator;

struct FinalizeCommunicatorEvent {
  FinalizeCommunicatorEvent(const StaticCommunicator & comm) : communicator(comm) {}
  const StaticCommunicator & communicator;
};

class CommunicatorEventHandler {
public:
  virtual ~CommunicatorEventHandler() {}
  virtual void onCommunicatorFinalize(__attribute__((unused)) const StaticCommunicator & communicator) { }
protected:
  inline void sendEvent(const FinalizeCommunicatorEvent & event) {
    onCommunicatorFinalize(event.communicator);
  }

  template<class EventHandler>
  friend class EventHandlerManager;
};


class StaticCommunicator : public EventHandlerManager<CommunicatorEventHandler>{
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
protected:
  StaticCommunicator(int & argc, char ** & argv,
		     CommunicatorType type = _communicator_mpi);

public:
  virtual ~StaticCommunicator() {
    FinalizeCommunicatorEvent *event = new FinalizeCommunicatorEvent(*this);
    this->sendEvent(*event);
 
   delete event;
    delete real_static_communicator;
  };

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Point to Point                                                           */
  /* ------------------------------------------------------------------------ */
  template<typename T> inline void send(T * buffer, Int size,
					Int receiver, Int tag);
  template<typename T> inline void receive(T * buffer, Int size,
					   Int sender, Int tag);

  template<typename T> inline CommunicationRequest * asyncSend(T * buffer,
							       Int size,
							       Int receiver,
							       Int tag);
  template<typename T> inline CommunicationRequest * asyncReceive(T * buffer,
								  Int size,
								  Int sender,
								  Int tag);

  template<typename T> inline void probe(Int sender, Int tag,
                                         CommunicationStatus & status);

  /* ------------------------------------------------------------------------ */
  /* Collectives                                                              */
  /* ------------------------------------------------------------------------ */
  template<typename T> inline void allReduce(T * values, Int nb_values,
					     const SynchronizerOperation & op);

  template<typename T> inline void allGather(T * values, Int nb_values);
  template<typename T> inline void allGatherV(T * values, Int * nb_values);

  template<typename T> inline void gather(T * values, Int nb_values,
					  Int root = 0);
  template<typename T> inline void gatherV(T * values, Int * nb_values,
					   Int root = 0);
  template<typename T> inline void broadcast(T * values, Int nb_values,
					     Int root = 0);

  inline void barrier();

  /* ------------------------------------------------------------------------ */
  /* Request handling                                                         */
  /* ------------------------------------------------------------------------ */
  inline bool testRequest(CommunicationRequest * request);

  inline void wait(CommunicationRequest * request);
  inline void waitAll(std::vector<CommunicationRequest *> & requests);

  inline void freeCommunicationRequest(CommunicationRequest * request);
  inline void freeCommunicationRequest(std::vector<CommunicationRequest *> & requests);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  virtual Int getNbProc() const { return real_static_communicator->psize; };
  virtual Int whoAmI() const { return real_static_communicator->prank; };

  AKANTU_GET_MACRO(RealStaticCommunicator, *real_static_communicator, const RealStaticCommunicator &);
  AKANTU_GET_MACRO_NOT_CONST(RealStaticCommunicator, *real_static_communicator, RealStaticCommunicator &);

  template<class Comm>
  Comm & getRealStaticCommunicator() { return dynamic_cast<Comm &>(*real_static_communicator); }

  static StaticCommunicator & getStaticCommunicator(CommunicatorType type = _communicator_mpi);

  static StaticCommunicator & getStaticCommunicator(int & argc, char ** & argv,
						    CommunicatorType type = _communicator_mpi);

  static bool isInstantiated() { return is_instantiated; };

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
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "static_communicator_inline_impl.hh"

/* -------------------------------------------------------------------------- */
/* Inline Functions ArrayBase                                                */
/* -------------------------------------------------------------------------- */
inline std::ostream & operator<<(std::ostream & stream, const CommunicationRequest & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__

#endif /* __AKANTU_STATIC_COMMUNICATOR_HH__ */
