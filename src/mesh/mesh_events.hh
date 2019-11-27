/**
 * @file   mesh_events.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Feb 20 2015
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Classes corresponding to mesh events type
 *
 * @section LICENSE
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include <utility>

#include "aka_array.hh"
#include "element.hh"
#include "element_type_map.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MESH_EVENTS_HH__
#define __AKANTU_MESH_EVENTS_HH__

namespace akantu {

/// akantu::MeshEvent is the base event for meshes
template <class Entity> class MeshEvent {
public:
  MeshEvent(const std::string & origin = "") : origin_(origin) {}

  virtual ~MeshEvent() = default;
  /// Get the list of entity modified by the event nodes or elements
  const Array<Entity> & getList() const { return list; }
  /// Get the list of entity modified by the event nodes or elements
  Array<Entity> & getList() { return list; }

  std::string origin() const { return origin_; }
  
protected:
  Array<Entity> list;

private:
  std::string origin_;
};

class Mesh;

/// akantu::MeshEvent related to new nodes in the mesh
class NewNodesEvent : public MeshEvent<UInt> {
public:
  NewNodesEvent(const std::string & origin = "") : MeshEvent(origin) {}
  ~NewNodesEvent() override = default;
};

/// akantu::MeshEvent related to nodes removed from the mesh
class RemovedNodesEvent : public MeshEvent<UInt> {
public:
  inline RemovedNodesEvent(const Mesh & mesh, const std::string & origin = "");

  ~RemovedNodesEvent() override = default;
  /// Get the new numbering following suppression of nodes from nodes arrays
  AKANTU_GET_MACRO_NOT_CONST(NewNumbering, new_numbering, Array<UInt> &);
  /// Get the new numbering following suppression of nodes from nodes arrays
  AKANTU_GET_MACRO(NewNumbering, new_numbering, const Array<UInt> &);

private:
  Array<UInt> new_numbering;
};

/// akantu::MeshEvent related to new elements in the mesh
class NewElementsEvent : public MeshEvent<Element> {
public:
  NewElementsEvent(const std::string & origin = "") : MeshEvent<Element>(origin) {}
  ~NewElementsEvent() override = default;
};

/// akantu::MeshEvent related to elements removed from the mesh
class RemovedElementsEvent : public MeshEvent<Element> {
public:
  inline RemovedElementsEvent(const Mesh & mesh,
                              const ID & new_numbering_id = "new_numbering",
                              const std::string & origin = "");

  ~RemovedElementsEvent() override = default;

  /// Get the new numbering following suppression of elements from elements
  /// arrays
  AKANTU_GET_MACRO(NewNumbering, new_numbering,
                   const ElementTypeMapArray<UInt> &);
  /// Get the new numbering following suppression of elements from elements
  /// arrays
  AKANTU_GET_MACRO_NOT_CONST(NewNumbering, new_numbering,
                             ElementTypeMapArray<UInt> &);
  /// Get the new numbering following suppression of elements from elements
  /// arrays
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(NewNumbering, new_numbering, UInt);
  /// Get the new numbering following suppression of elements from elements
  /// arrays
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(NewNumbering, new_numbering, UInt);

protected:
  ElementTypeMapArray<UInt> new_numbering;
};

/// akantu::MeshEvent for element that changed in some sort, can be seen as a
/// combination of removed and added elements
class ChangedElementsEvent : public RemovedElementsEvent {
public:
  inline ChangedElementsEvent(
      const Mesh & mesh, ID new_numbering_id = "changed_event:new_numbering",
      const std::string & origin = "")
      : RemovedElementsEvent(mesh, std::move(new_numbering_id), origin) {}

  ~ChangedElementsEvent() override = default;
  AKANTU_GET_MACRO(ListOld, list, const Array<Element> &);
  AKANTU_GET_MACRO_NOT_CONST(ListOld, list, Array<Element> &);
  AKANTU_GET_MACRO(ListNew, new_list, const Array<Element> &);
  AKANTU_GET_MACRO_NOT_CONST(ListNew, new_list, Array<Element> &);

protected:
  Array<Element> new_list;
};

/* -------------------------------------------------------------------------- */

class MeshEventHandler {
public:
  virtual ~MeshEventHandler() = default;
  /* ------------------------------------------------------------------------ */
  /* Internal code                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// send a akantu::NewNodesEvent
  inline void sendEvent(const NewNodesEvent & event) {
    onNodesAdded(event.getList(), event);
  }
  /// send a akantu::RemovedNodesEvent
  inline void sendEvent(const RemovedNodesEvent & event) {
    onNodesRemoved(event.getList(), event.getNewNumbering(), event);
  }
  /// send a akantu::NewElementsEvent
  inline void sendEvent(const NewElementsEvent & event) {
    onElementsAdded(event.getList(), event);
  }
  /// send a akantu::RemovedElementsEvent
  inline void sendEvent(const RemovedElementsEvent & event) {
    onElementsRemoved(event.getList(), event.getNewNumbering(), event);
  }
  /// send a akantu::ChangedElementsEvent
  inline void sendEvent(const ChangedElementsEvent & event) {
    onElementsChanged(event.getListOld(), event.getListNew(),
                      event.getNewNumbering(), event);
  }
  template <class EventHandler> friend class EventHandlerManager;

  /* ------------------------------------------------------------------------ */
  /* Interface                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// function to implement to react on  akantu::NewNodesEvent
  virtual void onNodesAdded(const Array<UInt> & /*nodes_list*/,
                            const NewNodesEvent & /*event*/) {}
  /// function to implement to react on  akantu::RemovedNodesEvent
  virtual void onNodesRemoved(const Array<UInt> & /*nodes_list*/,
                              const Array<UInt> & /*new_numbering*/,
                              const RemovedNodesEvent & /*event*/) {}
  /// function to implement to react on  akantu::NewElementsEvent
  virtual void onElementsAdded(const Array<Element> & /*elements_list*/,
                               const NewElementsEvent & /*event*/) {}
  /// function to implement to react on  akantu::RemovedElementsEvent
  virtual void
  onElementsRemoved(const Array<Element> & /*elements_list*/,
                    const ElementTypeMapArray<UInt> & /*new_numbering*/,
                    const RemovedElementsEvent & /*event*/) {}
  /// function to implement to react on  akantu::ChangedElementsEvent
  virtual void
  onElementsChanged(const Array<Element> & /*old_elements_list*/,
                    const Array<Element> & /*new_elements_list*/,
                    const ElementTypeMapArray<UInt> & /*new_numbering*/,
                    const ChangedElementsEvent & /*event*/) {}
};

} // namespace akantu

#endif /* __AKANTU_MESH_EVENTS_HH__ */
