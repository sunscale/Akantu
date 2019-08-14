/**
 * @file   communicator_event_handler.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Nov 15 2017
 *
 * @brief  Event handler of the communicator
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

#ifndef __AKANTU_COMMUNICATOR_EVENT_HANDLER_HH__
#define __AKANTU_COMMUNICATOR_EVENT_HANDLER_HH__

namespace akantu {

class Communicator;

struct FinalizeCommunicatorEvent {
  explicit FinalizeCommunicatorEvent(const Communicator & comm)
      : communicator(comm) {}
  const Communicator & communicator;
};

class CommunicatorEventHandler {
public:
  virtual ~CommunicatorEventHandler() = default;
  virtual void onCommunicatorFinalize() = 0;

private:
  inline void sendEvent(const FinalizeCommunicatorEvent &) {
    this->onCommunicatorFinalize();
  }

  template <class EventHandler> friend class EventHandlerManager;
};

} // namespace akantu

#endif /* __AKANTU_COMMUNICATOR_EVENT_HANDLER_HH__ */
