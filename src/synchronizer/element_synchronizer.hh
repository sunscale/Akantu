/**
 * @file   element_synchronizer.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Dana Christen <dana.christen@epfl.ch>
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Main element synchronizer
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

#ifndef AKANTU_ELEMENT_SYNCHRONIZER_HH_
#define AKANTU_ELEMENT_SYNCHRONIZER_HH_

/* -------------------------------------------------------------------------- */
#include "aka_array.hh"
#include "aka_common.hh"
#include "mesh_partition.hh"
#include "synchronizer_impl.hh"

namespace akantu {
class Mesh;
}

/* -------------------------------------------------------------------------- */
namespace akantu {

class ElementSynchronizer : public SynchronizerImpl<Element>,
                            public MeshEventHandler {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  ElementSynchronizer(Mesh & mesh, const ID & id = "element_synchronizer",
                      bool register_to_event_manager = true,
                      EventHandlerPriority event_priority = _ehp_synchronizer);

  ElementSynchronizer(const ElementSynchronizer & other,
                      const ID & id = "element_synchronizer_copy",
                      bool register_to_event_manager = true,
                      EventHandlerPriority event_priority = _ehp_synchronizer);

public:
  ~ElementSynchronizer() override;

  friend class ElementInfoPerProc;
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /// mesh event handler onElementsChanged
  void onElementsChanged(const Array<Element> & old_elements_list,
                         const Array<Element> & new_elements_list,
                         const ElementTypeMapArray<UInt> & new_numbering,
                         const ChangedElementsEvent & event) override;

  /// mesh event handler onRemovedElement
  void onElementsRemoved(const Array<Element> & element_to_remove,
                         const ElementTypeMapArray<UInt> & new_numbering,
                         const RemovedElementsEvent & event) override;

protected:
  /// remove elements from the synchronizer without renumbering them
  void removeElements(const Array<Element> & element_to_remove);

  /// renumber the elements in the synchronizer
  void renumberElements(const ElementTypeMapArray<UInt> & new_numbering);

  /// build processor to element correspondence
  void buildElementToPrank();

protected:
  /// fill the nodes type vector
  void fillNodesType(const MeshData & mesh_data,
                     DynamicCommunicationBuffer * buffers,
                     const std::string & tag_name, ElementType el_type,
                     const Array<UInt> & partition_num);

  template <typename T>
  void fillTagBufferTemplated(const MeshData & mesh_data,
                              DynamicCommunicationBuffer * buffers,
                              const std::string & tag_name,
                              ElementType el_type,
                              const Array<UInt> & partition_num,
                              const CSR<UInt> & ghost_partition);

  void fillTagBuffer(const MeshData & mesh_data,
                     DynamicCommunicationBuffer * buffers,
                     const std::string & tag_name, ElementType el_type,
                     const Array<UInt> & partition_num,
                     const CSR<UInt> & ghost_partition);

  /// function that handels the MeshData to be split (root side)
  static void synchronizeTagsSend(ElementSynchronizer & communicator, UInt root,
                                  Mesh & mesh, UInt nb_tags,
                                  ElementType type,
                                  const Array<UInt> & partition_num,
                                  const CSR<UInt> & ghost_partition,
                                  UInt nb_local_element, UInt nb_ghost_element);

  /// function that handles the MeshData to be split (other nodes)
  static void synchronizeTagsRecv(ElementSynchronizer & communicator, UInt root,
                                  Mesh & mesh, UInt nb_tags,
                                  ElementType type,
                                  UInt nb_local_element, UInt nb_ghost_element);

  /// function that handles the preexisting groups in the mesh
  static void synchronizeElementGroups(ElementSynchronizer & communicator,
                                       UInt root, Mesh & mesh,
                                       ElementType type,
                                       const Array<UInt> & partition_num,
                                       const CSR<UInt> & ghost_partition,
                                       UInt nb_element);

  /// function that handles the preexisting groups in the mesh
  static void synchronizeElementGroups(ElementSynchronizer & communicator,
                                       UInt root, Mesh & mesh,
                                       ElementType type);

  /// function that handles the preexisting groups in the mesh
  static void synchronizeNodeGroupsMaster(ElementSynchronizer & communicator,
                                          UInt root, Mesh & mesh);

  /// function that handles the preexisting groups in the mesh
  static void synchronizeNodeGroupsSlaves(ElementSynchronizer & communicator,
                                          UInt root, Mesh & mesh);

  template <class CommunicationBuffer>
  static void fillNodeGroupsFromBuffer(ElementSynchronizer & communicator,
                                       Mesh & mesh,
                                       CommunicationBuffer & buffer);

  /// substitute elements in the send and recv arrays
  void
  substituteElements(const std::map<Element, Element> & old_to_new_elements);

  /* ------------------------------------------------------------------------ */
  /* Sanity checks                                                            */
  /* ------------------------------------------------------------------------ */
  UInt sanityCheckDataSize(const Array<Element> & elements,
                           const SynchronizationTag & tag,
                           bool from_comm_desc = true) const override;
  void packSanityCheckData(CommunicationBuffer & /*buffer*/,
                           const Array<Element> & /*elements*/,
                           const SynchronizationTag & /*tag*/) const override;
  void unpackSanityCheckData(CommunicationBuffer & /*buffer*/,
                             const Array<Element> & /*elements*/,
                             const SynchronizationTag & /*tag*/, UInt /*proc*/,
                             UInt /*rank*/) const override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(Mesh, mesh, Mesh &);
  AKANTU_GET_MACRO(ElementToRank, element_to_prank,
                   const ElementTypeMapArray<Int> &);

  Int getRank(const Element & element) const final;
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// reference to the underlying mesh
  Mesh & mesh;

  friend class FacetSynchronizer;

  ElementTypeMapArray<Int> element_to_prank;
};

/* -------------------------------------------------------------------------- */
} // namespace akantu

#endif /* AKANTU_ELEMENT_SYNCHRONIZER_HH_ */
