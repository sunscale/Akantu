/**
 * @file   mesh_events.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Feb 20 10:48:36 2015
 *
 * @brief  Classes corresponding to mesh events type
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

#ifndef __AKANTU_MESH_EVENTS_HH__
#define __AKANTU_MESH_EVENTS_HH__

template<class Entity>
class MeshEvent {
public:
  virtual ~MeshEvent() {}
  const Array<Entity> & getList() const { return list; }
  Array<Entity> & getList() { return list; }
protected:
  Array<Entity> list;
};

class Mesh;

class NewNodesEvent : public MeshEvent<UInt> {
public:
  virtual ~NewNodesEvent() {};
};
class RemovedNodesEvent : public MeshEvent<UInt> {
public:
  virtual ~RemovedNodesEvent() {};
  inline RemovedNodesEvent(const Mesh & mesh);
  AKANTU_GET_MACRO_NOT_CONST(NewNumbering, new_numbering, Array<UInt> &);
  AKANTU_GET_MACRO(NewNumbering, new_numbering, const Array<UInt> &);
private:
  Array<UInt> new_numbering;
};

class NewElementsEvent : public MeshEvent<Element> {
public:
  virtual ~NewElementsEvent() {};
};

class RemovedElementsEvent : public MeshEvent<Element> {
public:
  virtual ~RemovedElementsEvent() {};
  inline RemovedElementsEvent(const Mesh & mesh);
  AKANTU_GET_MACRO(NewNumbering, new_numbering, const ElementTypeMapArray<UInt> &);
  AKANTU_GET_MACRO_NOT_CONST(NewNumbering, new_numbering, ElementTypeMapArray<UInt> &);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(NewNumbering, new_numbering, UInt);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(NewNumbering, new_numbering, UInt);
protected:
  ElementTypeMapArray<UInt> new_numbering;
};

class ChangedElementsEvent : public RemovedElementsEvent {
public:
  virtual ~ChangedElementsEvent() {};
  inline ChangedElementsEvent(const Mesh & mesh) : RemovedElementsEvent(mesh) {};
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
  virtual ~MeshEventHandler() {};
  /* ------------------------------------------------------------------------ */
  /* Internal code                                                            */
  /* ------------------------------------------------------------------------ */
private:
  inline void sendEvent(const NewNodesEvent & event)     { onNodesAdded  (event.getList(),
                                                                          event); }
  inline void sendEvent(const RemovedNodesEvent & event) { onNodesRemoved(event.getList(),
                                                                          event.getNewNumbering(),
                                                                          event); }

  inline void sendEvent(const NewElementsEvent & event)     { onElementsAdded  (event.getList(),
                                                                                event); }
  inline void sendEvent(const RemovedElementsEvent & event) { onElementsRemoved(event.getList(),
                                                                                event.getNewNumbering(),
                                                                                event); }

  inline void sendEvent(const ChangedElementsEvent & event) { onElementsChanged(event.getListOld(),
										event.getListNew(),
										event.getNewNumbering(),
										event); }

  template<class EventHandler>
  friend class EventHandlerManager;

  /* ------------------------------------------------------------------------ */
  /* Interface                                                                */
  /* ------------------------------------------------------------------------ */
public:
  virtual void onNodesAdded  (__attribute__((unused)) const Array<UInt> & nodes_list,
                              __attribute__((unused)) const NewNodesEvent & event) {  }
  virtual void onNodesRemoved(__attribute__((unused)) const Array<UInt> & nodes_list,
                              __attribute__((unused)) const Array<UInt> & new_numbering,
                              __attribute__((unused)) const RemovedNodesEvent & event) {  }

  virtual void onElementsAdded  (__attribute__((unused)) const Array<Element> & elements_list,
                                 __attribute__((unused)) const NewElementsEvent & event) { }
  virtual void onElementsRemoved(__attribute__((unused)) const Array<Element> & elements_list,
                                 __attribute__((unused)) const ElementTypeMapArray<UInt> & new_numbering,
                                 __attribute__((unused)) const RemovedElementsEvent & event) { }

  virtual void onElementsChanged(const Array<Element> & old_elements_list,
				 const Array<Element> & new_elements_list,
                                 const ElementTypeMapArray<UInt> & new_numbering,
                                 const ChangedElementsEvent & event) {
    NewElementsEvent new_elements;
    new_elements.getList() = new_elements_list;
    this->onElementsAdded(new_elements_list, new_elements);
    this->onElementsRemoved(old_elements_list, new_numbering, event);
  }
};


#endif /* __AKANTU_MESH_EVENTS_HH__ */
