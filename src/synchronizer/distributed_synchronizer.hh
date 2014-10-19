/**
 * @file   distributed_synchronizer.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Dana Christen <dana.christen@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Jun 16 2011
 * @date last modification: Fri Sep 05 2014
 *
 * @brief  wrapper to the static communicator
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

#ifndef __AKANTU_DISTRIBUTED_SYNCHRONIZER_HH__
#define __AKANTU_DISTRIBUTED_SYNCHRONIZER_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_array.hh"
#include "static_communicator.hh"
#include "synchronizer.hh"
#include "mesh.hh"
#include "mesh_partition.hh"
#include "communication_buffer.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class DistributedSynchronizer : public Synchronizer, public MeshEventHandler {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  DistributedSynchronizer(Mesh & mesh,
                          SynchronizerID id = "distributed_synchronizer",
                          MemoryID memory_id = 0);

public:
  virtual ~DistributedSynchronizer();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// get a  mesh and a partition and  create the local mesh  and the associated
  /// DistributedSynchronizer
  static DistributedSynchronizer *
  createDistributedSynchronizerMesh(Mesh & mesh,
                                    const MeshPartition * partition,
                                    UInt root = 0,
                                    SynchronizerID id = "distributed_synchronizer",
                                    MemoryID memory_id = 0);

  /* ------------------------------------------------------------------------ */
  /* Inherited from Synchronizer                                              */
  /* ------------------------------------------------------------------------ */

  /// asynchronous synchronization of ghosts
  void asynchronousSynchronize(DataAccessor & data_accessor,SynchronizationTag tag);

  /// wait end of asynchronous synchronization of ghosts
  void waitEndSynchronize(DataAccessor & data_accessor,SynchronizationTag tag);

  /// build processor to element corrispondance
  void buildPrankToElement();

  virtual void printself(std::ostream & stream, int indent = 0) const;

  /// mesh event handler onRemovedElement
  virtual void onElementsRemoved(const Array<Element> & element_list,
                                 const ElementTypeMapArray<UInt> & new_numbering,
                                 const RemovedElementsEvent & event);

  /// filter elements of a certain kind and copy them into a new synchronizer
  void filterElementsByKind(DistributedSynchronizer * new_synchronizer,
			    ElementKind kind);

  /// reset send and recv element lists
  void reset();

  /// compute buffer size for a given tag and data accessor
  void computeBufferSize(DataAccessor & data_accessor, SynchronizationTag tag);

protected:
  /// fill the nodes type vector
  void fillNodesType(Mesh & mesh);

  void fillNodesType(const MeshData & mesh_data,
                     DynamicCommunicationBuffer * buffers,
                     const std::string & tag_name,
                     const ElementType & el_type,
                     const Array<UInt> & partition_num);

  template<typename T>
  void fillTagBufferTemplated(const MeshData & mesh_data,
                              DynamicCommunicationBuffer * buffers,
                              const std::string & tag_name,
                              const ElementType & el_type,
                              const Array<UInt> & partition_num,
                              const CSR<UInt> & ghost_partition);

  void fillTagBuffer(const MeshData & mesh_data,
                     DynamicCommunicationBuffer * buffers,
                     const std::string & tag_name,
                     const ElementType & el_type,
                     const Array<UInt> & partition_num,
                     const CSR<UInt> & ghost_partition);

  template<typename T, typename BufferType>
  void populateMeshDataTemplated(MeshData & mesh_data,
                                 BufferType & buffer,
                                 const std::string & tag_name,
                                 const ElementType & el_type,
                                 UInt nb_component,
                                 UInt nb_local_element,
                                 UInt nb_ghost_element);

  template <typename BufferType>
  void populateMeshData(MeshData & mesh_data,
                        BufferType & buffer,
                        const std::string & tag_name,
                        const ElementType & el_type,
                        const MeshDataTypeCode & type_code,
                        UInt nb_component,
                        UInt nb_local_element,
                        UInt nb_ghost_element);

  /// fill the communications array of a distributedSynchronizer based on a partition array
  void fillCommunicationScheme(const UInt * partition,
                               UInt nb_local_element,
                               UInt nb_ghost_element,
                               ElementType type);

  /// function that handels the MeshData to be split (root side)
  static void synchronizeTagsSend(DistributedSynchronizer & communicator,
				  UInt root,
				  Mesh & mesh,
				  UInt nb_tags,
				  const ElementType & type,
				  const Array<UInt> & partition_num,
				  const CSR<UInt> & ghost_partition,
				  UInt nb_local_element,
				  UInt nb_ghost_element);

  /// function that handles the MeshData to be split (other nodes)
  static void synchronizeTagsRecv(DistributedSynchronizer & communicator,
				  UInt root,
				  Mesh & mesh,
				  UInt nb_tags,
				  const ElementType & type,
				  UInt nb_local_element,
				  UInt nb_ghost_element);

  /// function that handles the preexisting groups in the mesh
  static void synchronizeElementGroups(DistributedSynchronizer & communicator,
				       UInt root,
				       Mesh & mesh,
				       const ElementType & type,
				       const Array<UInt> & partition_num,
				       const CSR<UInt> & ghost_partition,
				       UInt nb_element);


  /// function that handles the preexisting groups in the mesh
  static void synchronizeElementGroups(DistributedSynchronizer & communicator,
				       UInt root,
				       Mesh & mesh,
				       const ElementType & type);

  template<class CommunicationBuffer>
  static void fillElementGroupsFromBuffer(DistributedSynchronizer & communicator,
					  Mesh & mesh,
					  const ElementType & type,
					  CommunicationBuffer & buffer);

  /// function that handles the preexisting groups in the mesh
  static void synchronizeNodeGroupsMaster(DistributedSynchronizer & communicator,
					  UInt root,
					  Mesh & mesh);


  /// function that handles the preexisting groups in the mesh
  static void synchronizeNodeGroupsSlaves(DistributedSynchronizer & communicator,
					  UInt root,
					  Mesh & mesh);

  template<class CommunicationBuffer>
  static void fillNodeGroupsFromBuffer(DistributedSynchronizer & communicator,
					  Mesh & mesh,
					  CommunicationBuffer & buffer);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  AKANTU_GET_MACRO(PrankToElement, prank_to_element, const ElementTypeMapArray<UInt> &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  enum CommTags {
    TAG_SIZES        = 0,
    TAG_CONNECTIVITY = 1,
    TAG_DATA         = 2,
    TAG_PARTITIONS   = 3,
    TAG_NB_NODES     = 4,
    TAG_NODES        = 5,
    TAG_COORDINATES  = 6,
    TAG_NODES_TYPE   = 7,
    TAG_MESH_DATA    = 8,
    TAG_ELEMENT_GROUP = 9,
    TAG_NODE_GROUP = 10,
  };

protected:
  /// reference to the underlying mesh
  Mesh & mesh;

  /// the static memory instance
  StaticCommunicator * static_communicator;

  class Communication {
  public:
    void resize(UInt size) {
      send_buffer.resize(size);
      recv_buffer.resize(size);
      size_to_send   .resize(size);
      size_to_receive.resize(size);
    }

  public:
    /// size of data to send to each processor
    std::vector<UInt> size_to_send;
    /// size of data to recv to each processor
    std::vector<UInt> size_to_receive;
    std::vector< CommunicationBuffer > send_buffer;
    std::vector< CommunicationBuffer > recv_buffer;

    std::vector<CommunicationRequest *> send_requests;
    std::vector<CommunicationRequest *> recv_requests;
  };

  std::map<SynchronizationTag, Communication> communications;

  /// list of element to send to proc p
  Array<Element> * send_element;
  /// list of element to receive from proc p
  Array<Element> * recv_element;

  UInt nb_proc;
  UInt rank;

  friend class FilteredSynchronizer;
  friend class FacetSynchronizer;

  ElementTypeMapArray<UInt> prank_to_element;

};

__END_AKANTU__

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "distributed_synchronizer_tmpl.hh"


#endif /* __AKANTU_DISTRIBUTED_SYNCHRONIZER_HH__ */
