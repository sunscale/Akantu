/**
 * @file   node_info_per_processor.hh
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
#include "mesh_accessor.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_NODE_INFO_PER_PROCESSOR_HH__
#define __AKANTU_NODE_INFO_PER_PROCESSOR_HH__

namespace akantu {
class NodeSynchronizer;
class StaticCommunicator;
}

/* -------------------------------------------------------------------------- */

namespace akantu {

class NodeInfoPerProc : protected MeshAccessor {
public:
  NodeInfoPerProc(NodeSynchronizer & synchronizer,
                  UInt message_cnt,
                  UInt root);

  virtual void synchronizeNodes() = 0;
  virtual void synchronizeTypes() = 0;
  virtual void synchronizeGroups() = 0;

protected:
  template <class CommunicationBuffer>
  void fillNodeGroupsFromBuffer(CommunicationBuffer & buffer);
  void fillNodesType();

  void fillCommunicationScheme(const Array<UInt> &);
protected:
  NodeSynchronizer & synchronizer;
  const StaticCommunicator & comm;
  UInt rank;
  UInt nb_proc;
  UInt root;

  Mesh & mesh;

  UInt spatial_dimension;

  UInt message_count;
};

/* -------------------------------------------------------------------------- */
class MasterNodeInfoPerProc : protected NodeInfoPerProc {
public:
  MasterNodeInfoPerProc(NodeSynchronizer & synchronizer,
                        UInt message_cnt,
                        UInt root);

  void synchronizeNodes();
  void synchronizeTypes();
  void synchronizeGroups();

private:
  /// get the list of nodes to send and send them
  std::vector<Array<UInt> > nodes_per_proc;
  Array<UInt> nb_nodes_per_proc;
};

/* -------------------------------------------------------------------------- */
class SlaveNodeInfoPerProc : protected NodeInfoPerProc {
public:
  SlaveNodeInfoPerProc(NodeSynchronizer & synchronizer,
                       UInt message_cnt,
                       UInt root);

  void synchronizeNodes();
  void synchronizeTypes();
  void synchronizeGroups();

private:
};

} // akantu

#endif /* __AKANTU_NODE_INFO_PER_PROCESSOR_HH__ */
