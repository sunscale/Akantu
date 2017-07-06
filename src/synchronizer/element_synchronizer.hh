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
 * @date last modification: Tue Dec 08 2015
 *
 * @brief  Main element synchronizer
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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

#ifndef __AKANTU_ELEMENT_SYNCHRONIZER_HH__
#define __AKANTU_ELEMENT_SYNCHRONIZER_HH__

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
  ElementSynchronizer(Mesh & mesh, const StaticCommunicator & comm,
                      const ID & id = "element_synchronizer",
                      MemoryID memory_id = 0,
                      bool register_to_event_manager = true,
                      UInt event_priority = 100);

  ElementSynchronizer(Mesh & mesh, const ID & id = "element_synchronizer",
                      MemoryID memory_id = 0,
                      bool register_to_event_manager = true,
                      UInt event_priority = 100);

public:
  virtual ~ElementSynchronizer();

  friend class ElementInfoPerProc;
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /// build processor to element correspondence
  void buildPrankToElement();

  /// mesh event handler onElementsChanged
  virtual void
  onElementsChanged(const Array<Element> & old_elements_list,
                    const Array<Element> & new_elements_list,
                    const ElementTypeMapArray<UInt> & new_numbering,
                    const ChangedElementsEvent & event);

  /// mesh event handler onRemovedElement
  virtual void
  onElementsRemoved(const Array<Element> & element_list,
                    const ElementTypeMapArray<UInt> & new_numbering,
                    const RemovedElementsEvent & event);
  /// mesh event handler onNodesAdded
  virtual void onNodesAdded(__attribute__((unused))
                            const Array<UInt> & nodes_list,
                            __attribute__((unused))
                            const NewNodesEvent & event){};

  /// mesh event handler onRemovedNodes
  virtual void
  onNodesRemoved(__attribute__((unused)) const Array<UInt> & nodes_list,
                 __attribute__((unused)) const Array<UInt> & new_numbering,
                 __attribute__((unused)) const RemovedNodesEvent & event){};

  /// mesh event handler onElementsAdded
  virtual void onElementsAdded(__attribute__((unused))
                               const Array<Element> & elements_list,
                               __attribute__((unused))
                               const NewElementsEvent & event){};

  /// filter elements of a certain kind and copy them into a new synchronizer
  void filterElementsByKind(ElementSynchronizer * new_synchronizer,
                            ElementKind kind);

  /// reset send and recv element lists
  void reset();

  /// remove elements from the synchronizer without renumbering them
  void removeElements(const Array<Element> & element_to_remove);

  /// renumber the elements in the synchronizer
  void renumberElements(const ElementTypeMapArray<UInt> & new_numbering);

protected:
  /// fill the nodes type vector
  void fillNodesType(const MeshData & mesh_data,
                     DynamicCommunicationBuffer * buffers,
                     const std::string & tag_name, const ElementType & el_type,
                     const Array<UInt> & partition_num);

  template <typename T>
  void fillTagBufferTemplated(const MeshData & mesh_data,
                              DynamicCommunicationBuffer * buffers,
                              const std::string & tag_name,
                              const ElementType & el_type,
                              const Array<UInt> & partition_num,
                              const CSR<UInt> & ghost_partition);

  void fillTagBuffer(const MeshData & mesh_data,
                     DynamicCommunicationBuffer * buffers,
                     const std::string & tag_name, const ElementType & el_type,
                     const Array<UInt> & partition_num,
                     const CSR<UInt> & ghost_partition);

  /// function that handels the MeshData to be split (root side)
  static void synchronizeTagsSend(ElementSynchronizer & communicator, UInt root,
                                  Mesh & mesh, UInt nb_tags,
                                  const ElementType & type,
                                  const Array<UInt> & partition_num,
                                  const CSR<UInt> & ghost_partition,
                                  UInt nb_local_element, UInt nb_ghost_element);

  /// function that handles the MeshData to be split (other nodes)
  static void synchronizeTagsRecv(ElementSynchronizer & communicator, UInt root,
                                  Mesh & mesh, UInt nb_tags,
                                  const ElementType & type,
                                  UInt nb_local_element, UInt nb_ghost_element);

  /// function that handles the preexisting groups in the mesh
  static void synchronizeElementGroups(ElementSynchronizer & communicator,
                                       UInt root, Mesh & mesh,
                                       const ElementType & type,
                                       const Array<UInt> & partition_num,
                                       const CSR<UInt> & ghost_partition,
                                       UInt nb_element);

  /// function that handles the preexisting groups in the mesh
  static void synchronizeElementGroups(ElementSynchronizer & communicator,
                                       UInt root, Mesh & mesh,
                                       const ElementType & type);

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
  virtual UInt sanityCheckDataSize(const Array<Element> & elements,
                                   const SynchronizationTag & tag) const;
  /* ------------------------------------------------------------------------ */
  virtual void
  packSanityCheckData(CommunicationDescriptor<Element> & comm_desc) const;
  /* ------------------------------------------------------------------------ */
  virtual void
  unpackSanityCheckData(CommunicationDescriptor<Element> & comm_desc) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(PrankToElement, prank_to_element,
                   const ElementTypeMapArray<UInt> &);

  AKANTU_GET_MACRO(Mesh, mesh, const Mesh &);
  AKANTU_GET_MACRO_NOT_CONST(Mesh, mesh, Mesh &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// reference to the underlying mesh
  Mesh & mesh;

  friend class FilteredSynchronizer;
  friend class FacetSynchronizer;

  ElementTypeMapArray<UInt> prank_to_element;
};

/* -------------------------------------------------------------------------- */
} // namespace akantu

#endif /* __AKANTU_ELEMENT_SYNCHRONIZER_HH__ */
