/**
 * @file   element_info_per_processor.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Mar 16 2016
 * @date last modification: Tue Nov 07 2017
 *
 * @brief  Helper classes to create the distributed synchronizer and distribute
 * a mesh
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
#include "aka_common.hh"
#include "communication_buffer.hh"
#include "mesh.hh"
#include "mesh_accessor.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_ELEMENT_INFO_PER_PROCESSOR_HH__
#define __AKANTU_ELEMENT_INFO_PER_PROCESSOR_HH__

namespace akantu {
class ElementSynchronizer;
class Communicator;
class MeshPartition;
} // namespace akantu

/* -------------------------------------------------------------------------- */
namespace akantu {

class ElementInfoPerProc : protected MeshAccessor {
public:
  ElementInfoPerProc(ElementSynchronizer & synchronizer, UInt message_cnt,
                     UInt root, ElementType type);
  bool synchronize();

protected:
  virtual void synchronizeConnectivities() = 0;
  virtual void synchronizePartitions() = 0;
  virtual void synchronizeTags() = 0;
  virtual void synchronizeGroups() = 0;
  virtual bool needSynchronize() = 0;

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

  UInt rank{0};
  UInt nb_proc{1};

  UInt root{0};

  ElementType type{_not_defined};

  UInt nb_tags{0};
  UInt nb_nodes_per_element{0};
  UInt nb_element{0};

  UInt nb_local_element{0};
  UInt nb_ghost_element{0};

  UInt message_count{0};
  Mesh & mesh;
  const Communicator & comm;
};

/* -------------------------------------------------------------------------- */
class MasterElementInfoPerProc : public ElementInfoPerProc {
public:
  MasterElementInfoPerProc(ElementSynchronizer & synchronizer, UInt message_cnt,
                           UInt root, ElementType type,
                           const MeshPartition & partition);

protected:
  void synchronizeConnectivities() override;
  void synchronizePartitions() override;
  void synchronizeTags() override;
  void synchronizeGroups() override;
  bool needSynchronize() override { return type != _not_defined; }

protected:
  template <typename T>
  void fillTagBufferTemplated(std::vector<DynamicCommunicationBuffer> & buffers,
                              const std::string & tag_name);
  void fillTagBuffer(std::vector<DynamicCommunicationBuffer> & buffers,
                     const std::string & tag_name);

private:
  const MeshPartition & partition;

  Vector<UInt> all_nb_local_element;
  Vector<UInt> all_nb_ghost_element;
  Vector<UInt> all_nb_element_to_send;
};

/* -------------------------------------------------------------------------- */
class SlaveElementInfoPerProc : public ElementInfoPerProc {
public:
  SlaveElementInfoPerProc(ElementSynchronizer & synchronizer, UInt message_cnt,
                          UInt root);

protected:
  void synchronizeConnectivities() override;
  void synchronizePartitions() override;
  void synchronizeTags() override;
  void synchronizeGroups() override;

  bool needSynchronize() override;

private:
  UInt nb_element_to_receive{0};
};

} // namespace akantu

#include "element_info_per_processor_tmpl.hh"

#endif /* __AKANTU_ELEMENT_INFO_PER_PROCESSOR_HH__ */
