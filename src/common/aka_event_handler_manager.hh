/**
 * @file   aka_event_handler_manager.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Sun Dec 03 2017
 *
 * @brief  Base of Event Handler classes
 *
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

#ifndef AKANTU_AKA_EVENT_HANDLER_MANAGER_HH_
#define AKANTU_AKA_EVENT_HANDLER_MANAGER_HH_

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
#include <list>
/* -------------------------------------------------------------------------- */

namespace akantu {

template <class EventHandler> class EventHandlerManager {
private:
  using priority_value = std::pair<EventHandlerPriority, EventHandler *>;
  using priority_list = std::list<priority_value>;
  struct KeyComp {
    bool operator()(const priority_value & a, const priority_value & b) const {
      return (a.first < b.first);
    }
    bool operator()(const priority_value & a, UInt b) const {
      return (a.first < b);
    }
  };

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  virtual ~EventHandlerManager() = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// register a new EventHandler to the Manager. The register object
  /// will then be informed about the events the manager observes.
  void registerEventHandler(EventHandler & event_handler,
                            EventHandlerPriority priority = _ehp_highest) {
    auto it = this->searchEventHandler(event_handler);

    if (it != this->event_handlers.end()) {
      AKANTU_EXCEPTION("This event handler was already registered (priority: "
                       << priority << ")");
    }

    auto pos =
        std::lower_bound(this->event_handlers.begin(),
                         this->event_handlers.end(), priority, KeyComp());

    this->event_handlers.insert(pos, std::make_pair(priority, &event_handler));
  }

  /// unregister a EventHandler object. This object will not be
  /// notified anymore about the events this manager observes.
  void unregisterEventHandler(EventHandler & event_handler) {
    auto it = this->searchEventHandler(event_handler);

    if (it == this->event_handlers.end()) {
      AKANTU_EXCEPTION("This event handler is not registered");
    }

    this->event_handlers.erase(it);
  }

  /// Notify all the registered EventHandlers about the event that just occured.
  template <class Event> void sendEvent(const Event & event) {
    for (auto & pair : this->event_handlers) {
      pair.second->sendEvent(event);
    }
  }

private:
  typename priority_list::iterator searchEventHandler(EventHandler & handler) {
    auto it = this->event_handlers.begin();
    auto end = this->event_handlers.end();

    for (; it != end && it->second != &handler; ++it) {
      ;
    }

    return it;
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// list of the event handlers
  priority_list event_handlers;
};

} // namespace akantu

#endif /* AKANTU_AKA_EVENT_HANDLER_MANAGER_HH_ */
