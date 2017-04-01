/**
 * @file   element_info_per_processor.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Mar 11 14:45:15 2016
 *
 * @brief  Helper classes to create the distributed synchronizer and distribute
 *         a mesh
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
#include "aka_common.hh"
#include "communication_buffer.hh"
#include "mesh.hh"
#include "mesh_accessor.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_ELEMENT_INFO_PER_PROCESSOR_HH__
#define __AKANTU_ELEMENT_INFO_PER_PROCESSOR_HH__

namespace akantu {
class ElementSynchronizer;
class StaticCommunicator;
class MeshPartition;
}

/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__

class ElementInfoPerProc : protected MeshAccessor {
public:
  ElementInfoPerProc(ElementSynchronizer & synchronizer, UInt message_cnt,
                     UInt root, ElementType type);

  virtual void synchronizeConnectivities() = 0;
  virtual void synchronizePartitions() = 0;
  virtual void synchronizeTags() = 0;
  virtual void synchronizeGroups() = 0;

protected:
  void fillCommunicationScheme(const Array<UInt> & partition);

  template <class CommunicationBuffer>
  void fillElementGroupsFromBuffer(CommunicationBuffer & buffer);

  template <typename T, typename BufferType>
  void fillMeshDataTemplated(BufferType & buffer, const std::string & tag_name,
                             UInt nb_component);

  template <typename BufferType>
  void fillMeshData(BufferType & buffer, const std::string & tag_name,
                    const MeshDataTypeCode & type_code, UInt nb_component);

protected:
  ElementSynchronizer & synchronizer;

  UInt rank;
  UInt nb_proc;

  UInt root;

  ElementType type;

  UInt nb_tags;
  UInt nb_nodes_per_element;
  UInt nb_element;

  UInt nb_local_element;
  UInt nb_ghost_element;

  UInt message_count;
  Mesh & mesh;
  const StaticCommunicator & comm;
};

/* -------------------------------------------------------------------------- */
class MasterElementInfoPerProc : protected ElementInfoPerProc {
public:
  MasterElementInfoPerProc(ElementSynchronizer & synchronizer, UInt message_cnt,
                           UInt root, ElementType type,
                           const MeshPartition & partition);

  void synchronizeConnectivities();
  void synchronizePartitions();
  void synchronizeTags();
  void synchronizeGroups();

protected:
  template <typename T>
  void fillTagBufferTemplated(DynamicCommunicationBuffer * buffers,
                              const std::string & tag_name);
  void fillTagBuffer(DynamicCommunicationBuffer * buffers,
                     const std::string & tag_name);

private:
  const MeshPartition & partition;

  Vector<UInt> all_nb_local_element;
  Vector<UInt> all_nb_ghost_element;
  Vector<UInt> all_nb_element_to_send;
};

/* -------------------------------------------------------------------------- */
class SlaveElementInfoPerProc : protected ElementInfoPerProc {
public:
  SlaveElementInfoPerProc(ElementSynchronizer & synchronizer, UInt message_cnt,
                          UInt root);

  void synchronizeConnectivities();
  void synchronizePartitions();
  void synchronizeTags();
  void synchronizeGroups();

  bool needSynchronize();

private:
  UInt nb_element_to_receive;
};

__END_AKANTU__

#include "element_info_per_processor_tmpl.hh"

#endif /* __AKANTU_ELEMENT_INFO_PER_PROCESSOR_HH__ */
