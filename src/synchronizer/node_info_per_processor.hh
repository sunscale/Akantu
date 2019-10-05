/**
 * @file   node_info_per_processor.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Mar 16 2016
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  Helper classes to create the distributed synchronizer and distribute
 * a mesh
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
#include "communication_buffer.hh"
#include "mesh_accessor.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_NODE_INFO_PER_PROCESSOR_HH__
#define __AKANTU_NODE_INFO_PER_PROCESSOR_HH__

namespace akantu {
class NodeSynchronizer;
class Communicator;
} // namespace akantu

/* -------------------------------------------------------------------------- */

namespace akantu {

class NodeInfoPerProc : protected MeshAccessor {
public:
  NodeInfoPerProc(NodeSynchronizer & synchronizer, UInt message_cnt, UInt root);

  void synchronize();

protected:
  virtual void synchronizeNodes() = 0;
  virtual void synchronizeTypes() = 0;
  virtual void synchronizeGroups() = 0;
  virtual void synchronizePeriodicity() = 0;
  virtual void synchronizeTags() = 0;

protected:
  template <class CommunicationBuffer>
  void fillNodeGroupsFromBuffer(CommunicationBuffer & buffer);
  void fillNodesType();

  void fillCommunicationScheme(const Array<UInt> &);
  void fillNodalData(DynamicCommunicationBuffer & buffer, std::string tag_name);

  void fillPeriodicPairs(const Array<UInt> &, std::vector<UInt> &);
  void receiveMissingPeriodic(DynamicCommunicationBuffer &);

protected:
  NodeSynchronizer & synchronizer;
  const Communicator & comm;
  UInt rank;
  UInt nb_proc;
  UInt root;

  Mesh & mesh;

  UInt spatial_dimension;
  UInt message_count;
};

/* -------------------------------------------------------------------------- */
class MasterNodeInfoPerProc : public NodeInfoPerProc {
public:
  MasterNodeInfoPerProc(NodeSynchronizer & synchronizer, UInt message_cnt,
                        UInt root);

  void synchronizeNodes() override;
  void synchronizeTypes() override;
  void synchronizeGroups() override;
  void synchronizePeriodicity() override;
  void synchronizeTags() override;

private:
  void fillTagBuffers(std::vector<DynamicCommunicationBuffer> & buffers,
                      const std::string & tag_name);

  /// get the list of nodes to send and send them
  std::vector<Array<UInt>> nodes_per_proc;
  Array<UInt> nb_nodes_per_proc;
  Array<Real> all_nodes;
  Array<NodeFlag> all_periodic_flags;
  Array<Int> nodes_pranks;
};

/* -------------------------------------------------------------------------- */
class SlaveNodeInfoPerProc : public NodeInfoPerProc {
public:
  SlaveNodeInfoPerProc(NodeSynchronizer & synchronizer, UInt message_cnt,
                       UInt root);

  void synchronizeNodes() override;
  void synchronizeTypes() override;
  void synchronizeGroups() override;
  void synchronizePeriodicity() override;
  void synchronizeTags() override;

private:
};

} // namespace akantu

#endif /* __AKANTU_NODE_INFO_PER_PROCESSOR_HH__ */
