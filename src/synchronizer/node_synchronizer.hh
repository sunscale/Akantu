/**
 * @file   node_synchronizer.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Nov 08 2016
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Synchronizer for nodal information
 *
 *
 * Copyright (©) 2016-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "mesh_events.hh"
#include "synchronizer_impl.hh"
/* -------------------------------------------------------------------------- */
#include <unordered_map>
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

  ~NodeSynchronizer() override;

  UInt sanityCheckDataSize(const Array<UInt> & nodes,
                           const SynchronizationTag & tag,
                           bool from_comm_desc) const override;
  void packSanityCheckData(CommunicationBuffer & buffer,
                           const Array<UInt> & nodes,
                           const SynchronizationTag & /*tag*/) const override;
  void unpackSanityCheckData(CommunicationBuffer & buffer,
                             const Array<UInt> & nodes,
                             const SynchronizationTag & tag, UInt proc,
                             UInt rank) const override;

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

  /* ------------------------------------------------------------------------ */
  NodeSynchronizer & operator=(const NodeSynchronizer & other) {
    copySchemes(other);
    return *this;
  }

  friend class NodeInfoPerProc;
protected:
  void fillEntityToSend(Array<UInt> & entities_to_send) override;

public:
  AKANTU_GET_MACRO(Mesh, mesh, Mesh &);

  inline UInt canScatterSize() override;
  inline UInt gatheredSize() override;

  inline UInt localToGlobalEntity(const UInt & local) override;
  
protected:
  Int getRank(const UInt & node) const final;

protected:
  Mesh & mesh;
};

} // namespace akantu

#include "node_synchronizer_inline_impl.hh"

#endif /* __AKANTU_NODE_SYNCHRONIZER_HH__ */
