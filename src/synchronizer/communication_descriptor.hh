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
#include "communication_tag.hh"
#include "data_accessor.hh"
#include "communication_request.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_COMMUNICATION_DESCRIPTOR_HH__
#define __AKANTU_COMMUNICATION_DESCRIPTOR_HH__

namespace akantu {

/* ------------------------------------------------------------------------ */
enum CommunicationSendRecv  { _send, _recv, _csr_not_defined };

/* -------------------------------------------------------------------------- */
struct CommunicationSRType {
  typedef CommunicationSendRecv type;
  static const type _begin_ = _send;
  static const type _end_ = _csr_not_defined;
};

typedef safe_enum<CommunicationSRType> send_recv_t;

namespace {
  send_recv_t iterate_send_recv{};
}

/* ------------------------------------------------------------------------ */
class Communication {
public:
  explicit Communication(const CommunicationSendRecv & type = _csr_not_defined)
      : _size(0), _type(type) {}

  Communication(const Communication &) = delete;
  Communication & operator=(const Communication &) = delete;

  void resize(UInt size) {
    this->_size = size;
    this->_buffer.resize(size);
  }

  inline const CommunicationSendRecv & type() const { return this->_type; }
  inline const UInt & size() const { return this->_size; }

  inline const CommunicationRequest & request() const { return this->_request; }
  inline CommunicationRequest & request() { return this->_request; }

  inline const CommunicationBuffer & buffer() const { return this->_buffer; }
  inline CommunicationBuffer & buffer() { return this->_buffer; }

private:
  UInt _size;
  CommunicationBuffer _buffer;
  CommunicationRequest _request;
  CommunicationSendRecv _type;
};

template <class Entity> class Communications;

/* ------------------------------------------------------------------------ */
template <class Entity> class CommunicationDescriptor {
public:
  CommunicationDescriptor(Communication & communication, Array<Entity> & scheme,
                          Communications<Entity> & communications,
                          const SynchronizationTag & tag, UInt proc);

  CommunicationDescriptor(const CommunicationDescriptor &) = default;

  CommunicationDescriptor &
  operator=(const CommunicationDescriptor &) = default;

  /// get the quantity of data in the buffer
  UInt getNbData() { return communication.size(); }
  /// set the quantity of data in the buffer
  void setNbData(UInt size) { communication.resize(size); }

  /// get the corresponding tag
  const SynchronizationTag & getTag() const { return tag; }
  /// get the data buffer
  CommunicationBuffer & getBuffer();

  /// get the corresponding request
  CommunicationRequest & getRequest();

  /// get the communication scheme
  Array<Entity> & getScheme();

  /// reset the buffer before pack or after unpack
  void resetBuffer();

  /// pack data for entities in the buffer
  void packData(const DataAccessor<Entity> & accessor);

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
} // namespace akantu

#include "communication_descriptor_tmpl.hh"

#endif /* __AKANTU_COMMUNICATION_DESCRIPTOR_HH__ */
