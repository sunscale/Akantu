/**
 * @file   test_node_synchronizer.cc
 *
 * @author Nicolas Richart
 *
 * @date creation  Tue Apr 04 2017
 *
 * @brief test the default node synchronizer present in the mesh
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
#include "data_accessor.hh"
#include "mesh.hh"
#include "node_synchronizer.hh"
#include "communicator.hh"
/* -------------------------------------------------------------------------- */
#include <limits>
/* -------------------------------------------------------------------------- */

using namespace akantu;

class DataAccessorTest : public DataAccessor<UInt> {
public:
  explicit DataAccessorTest(Array<int> & data) : data(data) {}

  UInt getNbData(const Array<UInt> & nodes,
                 const SynchronizationTag &) const {
    return nodes.size() * sizeof(int);
  }

  void packData(CommunicationBuffer & buffer, const Array<UInt> & nodes,
                const SynchronizationTag &) const {
    for (auto node : nodes) {
      buffer << data(node);
    }
  }

  void unpackData(CommunicationBuffer & buffer, const Array<UInt> & nodes,
                  const SynchronizationTag &) {
    for (auto node : nodes) {
      buffer >> data(node);
    }
  }

protected:
  Array<int> & data;
};

int main(int argc, char * argv[]) {
  initialize(argc, argv);

  UInt spatial_dimension = 3;

  Mesh mesh(spatial_dimension);

  const auto & comm = Communicator::getStaticCommunicator();
  Int prank = comm.whoAmI();

  if (prank == 0)
    mesh.read("cube.msh");

  mesh.distribute();

  UInt nb_nodes = mesh.getNbNodes();
  Array<int> node_data(nb_nodes);

  constexpr int max_int = std::numeric_limits<int>::max();

  for (UInt n = 0; n < nb_nodes; ++n) {
    UInt gn = mesh.getNodeGlobalId(n);
    if (mesh.isMasterNode(n))
      node_data(n) = gn;
    else if (mesh.isLocalNode(n))
      node_data(n) = -gn;
    else if (mesh.isSlaveNode(n))
      node_data(n) = max_int;
    else
      node_data(n) = -max_int;
  }

  DataAccessorTest data_accessor(node_data);

  auto & node_synchronizer = mesh.getNodeSynchronizer();
  node_synchronizer.synchronize(data_accessor, _gst_test);


  for (UInt n = 0; n < nb_nodes; ++n) {
    int gn = mesh.getNodeGlobalId(n);
    if(!((mesh.isMasterNode(n) && node_data(n) == gn) ||
         (mesh.isLocalNode(n)  && node_data(n) == -gn) ||
         (mesh.isSlaveNode(n)  && node_data(n) == gn) ||
         (mesh.isPureGhostNode(n) && node_data(n) == -max_int)
         )
      )
      {
        debug::setDebugLevel(dblTest);
        std::cout << "prank : " << prank << " (node " << gn << "[" << n << "]) - "
                  << "( isMaster: " << mesh.isMasterNode(n) << " && " << node_data(n) << " == "  << gn << ") "
                  << "( isLocal: " << mesh.isLocalNode(n) << " && " << node_data(n) << " == "  <<  - gn << ") "
                  << "( isSlave: " << mesh.isSlaveNode(n) << " && " << node_data(n) << " == "  << gn << ") "
                  << "( isPureGhost: " << mesh.isPureGhostNode(n) << " && " << node_data(n) << " == "  << -max_int << ") "
                  << std::endl;
        debug::setDebugLevel(dblDebug);
        return -1;
      }
  }

  finalize();

  return 0;
}
