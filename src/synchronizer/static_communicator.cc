/**
 * @file   static_communicator.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Jan 19 2016
 *
 * @brief  implementation of the common part of the static communicator
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
#include "static_communicator.hh"
#include "static_communicator_dummy.hh"

#ifdef AKANTU_USE_MPI
#include "static_communicator_mpi.hh"
#endif

/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
bool StaticCommunicator::is_instantiated = false;
StaticCommunicator * StaticCommunicator::static_communicator = NULL;

UInt InternalCommunicationRequest::counter = 0;

/* -------------------------------------------------------------------------- */
InternalCommunicationRequest::InternalCommunicationRequest(UInt source, UInt dest)
    : source(source), destination(dest) {
  this->id = counter++;
}

/* -------------------------------------------------------------------------- */
InternalCommunicationRequest::~InternalCommunicationRequest() {}

/* -------------------------------------------------------------------------- */
void InternalCommunicationRequest::printself(std::ostream & stream, int indent) const {
  std::string space;
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
    ;

  stream << space << "CommunicationRequest [" << std::endl;
  stream << space << " + id          : " << id << std::endl;
  stream << space << " + source      : " << source << std::endl;
  stream << space << " + destination : " << destination << std::endl;
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
StaticCommunicator::StaticCommunicator(int & argc, char **& argv,
                                       CommunicatorType type) {
  real_type = type;
#ifdef AKANTU_USE_MPI
  if (type == _communicator_mpi) {
    real_static_communicator = new StaticCommunicatorMPI(argc, argv);
  } else {
#endif
    real_static_communicator = new StaticCommunicatorDummy(argc, argv);
#ifdef AKANTU_USE_MPI
  }
#endif
  this->is_instantiated = true;
}

/* -------------------------------------------------------------------------- */
StaticCommunicator::~StaticCommunicator() {
  FinalizeCommunicatorEvent * event = new FinalizeCommunicatorEvent(*this);
  this->sendEvent(*event);

  this->barrier();
  delete event;
  delete real_static_communicator;
  is_instantiated = false;
  StaticCommunicator::static_communicator = NULL;
}

/* -------------------------------------------------------------------------- */
StaticCommunicator &
StaticCommunicator::getStaticCommunicator(CommunicatorType type) {
  AKANTU_DEBUG_IN();

  if (!static_communicator) {
    int nb_args = 0;
    char ** null;
    static_communicator = new StaticCommunicator(nb_args, null, type);
  }


  AKANTU_DEBUG_OUT();
  return *static_communicator;
}

/* -------------------------------------------------------------------------- */
StaticCommunicator &
StaticCommunicator::getStaticCommunicator(int & argc, char **& argv,
                                          CommunicatorType type) {
  if (!static_communicator)
    static_communicator = new StaticCommunicator(argc, argv, type);

  return getStaticCommunicator(type);
}

/* -------------------------------------------------------------------------- */
StaticCommunicator * StaticCommunicator::getStaticCommunicatorDummy() {
  int nb_args = 0;
  char ** null;
  return new StaticCommunicator(nb_args, null, _communicator_dummy);
}


__END_AKANTU__
