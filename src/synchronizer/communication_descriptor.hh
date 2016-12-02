/**
 * @file   synchronizer_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Aug  3 13:49:36 2016
 *
 * @brief  Implementation of the helper classes for the synchronizer
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_array.hh"
#include "data_accessor.hh"
#include "communication_tag.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_COMMUNICATION_DESCRIPTOR_HH__
#define __AKANTU_COMMUNICATION_DESCRIPTOR_HH__

namespace akantu {

/* ------------------------------------------------------------------------ */
enum CommunicationSendRecv {
  _send,
  _recv,
  _csr_not_defined
};

/* ------------------------------------------------------------------------ */
struct Communication {
  Communication() : size(0), type(_csr_not_defined) {}
  UInt size;
  CommunicationBuffer buffer;
  CommunicationRequest request;
  CommunicationSendRecv type;
};

template <class Entity> class Communications;

/* ------------------------------------------------------------------------ */
template <class Entity> class CommunicationDescriptor {
public:
  CommunicationDescriptor(Communication & communication, Array<Entity> & scheme,
                          Communications<Entity> & communications,
                          const SynchronizationTag & tag, UInt proc);
  /// get the quantity of data in the buffer
  UInt getNbData() { return communication.size; }
  /// set the quantity of data in the buffer
  void setNbData(UInt size) { communication.size = size; }

  /// get the corresponding tag
  const SynchronizationTag & getTag() const { return tag; }
  /// get the data buffer
  CommunicationBuffer & getBuffer();

  /// get the corresponding request
  CommunicationRequest & getRequest();

  /// get the communication scheme
  Array<Entity> & getScheme();

  /// pack data for entities in the buffer
  void packData(DataAccessor<Entity> & accessor);

  /// unpack data for entities from the buffer
  void unpackData(DataAccessor<Entity> & accessor);

  /// posts asynchronous send requests
  void postSend(int hash_id);

  /// posts asynchronous receive requests
  void postRecv(int hash_id);

  /// free the request
  void freeRequest();

  UInt getProc() { return proc; }

protected:
  Communication & communication;
  Array<Entity> & scheme;
  Communications<Entity> & communications;
  const SynchronizationTag & tag;
  UInt proc;
  UInt rank;
  UInt counter;
};

/* -------------------------------------------------------------------------- */
} // akantu

#include "communication_descriptor_tmpl.hh"

#endif /* __AKANTU_COMMUNICATION_DESCRIPTOR_HH__ */
