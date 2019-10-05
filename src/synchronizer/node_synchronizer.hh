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
 * @section LICENSE
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

  /* ------------------------------------------------------------------------ */
  /// Uses the synchronizer to perform a reduction on the vector
  template <template <class> class Op, typename T>
  void reduceSynchronize(Array<T> & array) const;

  /* ------------------------------------------------------------------------ */
  template <typename T> void synchronizeData(Array<T> & array) const;

  friend class NodeInfoPerProc;

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

public:
  AKANTU_GET_MACRO(Mesh, mesh, Mesh &);

protected:
  Int getRank(const UInt & node) const final;

protected:
  Mesh & mesh;

  // std::unordered_map<UInt, Int> node_to_prank;
};

/* -------------------------------------------------------------------------- */
template <typename T>
void NodeSynchronizer::synchronizeData(Array<T> & array) const {
  SimpleUIntDataAccessor<T> data_accessor(array, SynchronizationTag::_whatever);
  this->synchronizeOnce(data_accessor, SynchronizationTag::_whatever);
}

/* -------------------------------------------------------------------------- */
template <template <class> class Op, typename T>
void NodeSynchronizer::reduceSynchronize(Array<T> & array) const {
  ReduceDataAccessor<UInt, Op, T> data_accessor(array,
                                                SynchronizationTag::_whatever);
  this->slaveReductionOnceImpl(data_accessor, SynchronizationTag::_whatever);
  this->synchronizeData(array);
}

} // namespace akantu

#endif /* __AKANTU_NODE_SYNCHRONIZER_HH__ */
