/**
 * @file   aka_event_handler_manager.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Aug 23 2012
 * @date last modification: Mon Jun 02 2014
 *
 * @brief  Base of Event Handler classes
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

#ifndef __AKANTU_AKA_EVENT_HANDLER_MANAGER_HH__
#define __AKANTU_AKA_EVENT_HANDLER_MANAGER_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */
#include <set>
/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

template<class EventHandler>
class EventHandlerManager {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  virtual ~EventHandlerManager() {};
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  void registerEventHandler(EventHandler & event_handler) {
    event_handlers.insert(&event_handler);
  }

  void unregisterEventHandler(EventHandler & event_handler) {
    event_handlers.erase(&event_handler);
  }

  template<class Event>
  void sendEvent(const Event & event) {
    typename std::set<EventHandler *>::iterator it = event_handlers.begin();
    typename std::set<EventHandler *>::iterator end = event_handlers.end();
    for(;it != end; ++it)
      (*it)->sendEvent(event);
  }

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// list of the event handlers
  std::set<EventHandler *> event_handlers;
};

__END_AKANTU__

#endif /* __AKANTU_AKA_EVENT_HANDLER_MANAGER_HH__ */
