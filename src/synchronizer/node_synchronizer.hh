/**
 * @file   node_synchronizer.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Sep 23 11:45:24 2016
 *
 * @brief  Synchronizer for nodal information
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
#include "mesh_events.hh"
#include "synchronizer_impl.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_NODE_SYNCHRONIZER_HH__
#define __AKANTU_NODE_SYNCHRONIZER_HH__

namespace akantu {

class NodeSynchronizer : public MeshEventHandler,
                         public SynchronizerImpl<UInt> {
public:
  NodeSynchronizer(Mesh & mesh, const ID & id = "element_synchronizer",
                   MemoryID memory_id = 0,
                   const bool register_to_event_manager = true,
                   EventHandlerPriority event_priority = _ehp_synchronizer);

  virtual ~NodeSynchronizer();

  friend class NodeInfoPerProc;

  /// function to implement to react on  akantu::NewNodesEvent
  void onNodesAdded(const Array<UInt> &, const NewNodesEvent &) override;

  /// function to implement to react on  akantu::RemovedNodesEvent
  void onNodesRemoved(const Array<UInt> &, const Array<UInt> &,
                      const RemovedNodesEvent &) override {}
  /// function to implement to react on  akantu::NewElementsEvent
  void onElementsAdded(const Array<Element> &,
                       const NewElementsEvent &) override {}
  /// function to implement to react on  akantu::RemovedElementsEvent
  void onElementsRemoved(const Array<Element> &,
                         const ElementTypeMapArray<UInt> &,
                         const RemovedElementsEvent &) override {}
  /// function to implement to react on  akantu::ChangedElementsEvent
  void onElementsChanged(const Array<Element> &, const Array<Element> &,
                         const ElementTypeMapArray<UInt> &,
                         const ChangedElementsEvent &) override {}

public:
  AKANTU_GET_MACRO(Mesh, mesh, Mesh &);

protected:
  Mesh & mesh;
};

} // namespace akantu

#endif /* __AKANTU_NODE_SYNCHRONIZER_HH__ */
