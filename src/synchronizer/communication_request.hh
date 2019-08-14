/**
 * @file   communication_request.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Jun 14 2010
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  empty class just for inheritance
 *
 * @section LICENSE
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
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */
#include <memory>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_REAL_STATIC_COMMUNICATOR_HH__
#define __AKANTU_REAL_STATIC_COMMUNICATOR_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
class InternalCommunicationRequest {
public:
  InternalCommunicationRequest(UInt source, UInt dest);
  virtual ~InternalCommunicationRequest();

  virtual void printself(std::ostream & stream, int indent = 0) const;

  AKANTU_GET_MACRO(Source, source, UInt);
  AKANTU_GET_MACRO(Destination, destination, UInt);

private:
  UInt source;
  UInt destination;
  UInt id;
  static UInt counter;
};

/* -------------------------------------------------------------------------- */
class CommunicationRequest {
public:
  CommunicationRequest(
      std::shared_ptr<InternalCommunicationRequest> request = nullptr)
      : request(std::move(request)) {}

  virtual ~CommunicationRequest() = default;

  virtual void free() { request.reset(); }

  void printself(std::ostream & stream, int indent = 0) const {
    request->printself(stream, indent);
  };

  UInt getSource() const { return request->getSource(); }
  UInt getDestination() const { return request->getDestination(); }

  bool isFreed() const { return request.get() == nullptr; }

  InternalCommunicationRequest & getInternal() { return *request; }

private:
  std::shared_ptr<InternalCommunicationRequest> request;
};

/* -------------------------------------------------------------------------- */
class CommunicationStatus {
public:
  AKANTU_GET_MACRO(Source, source, Int);
  UInt size() const { return size_; }
  AKANTU_GET_MACRO(Tag, tag, Int);

  AKANTU_SET_MACRO(Source, source, Int);
  AKANTU_SET_MACRO(Size, size_, UInt);
  AKANTU_SET_MACRO(Tag, tag, Int);

private:
  Int source{0};
  UInt size_{0};
  Int tag{0};
};

/* -------------------------------------------------------------------------- */
/// Datatype to pack pairs for MPI_{MIN,MAX}LOC
template <typename T1, typename T2> struct SCMinMaxLoc {
  T1 min_max;
  T2 loc;
};

} // namespace akantu

#endif /* __AKANTU_REAL_STATIC_COMMUNICATOR_HH__ */
