/**
 * @file   test_node_synchronizer.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu May 11 2017
 * @date last modification: Fri Jan 26 2018
 *
 * @brief  test the default node synchronizer present in the mesh
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
#include "test_synchronizers_fixture.hh"
/* -------------------------------------------------------------------------- */
#include "communicator.hh"
#include "data_accessor.hh"
#include "mesh.hh"
#include "node_synchronizer.hh"
/* -------------------------------------------------------------------------- */
#include <chrono>
#include <limits>
#include <random>
#include <thread>
/* -------------------------------------------------------------------------- */

using namespace akantu;

class DataAccessorTest : public DataAccessor<UInt> {
public:
  explicit DataAccessorTest(Array<int> & data) : data(data) {}

  UInt getNbData(const Array<UInt> & nodes, const SynchronizationTag &) const {
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

class TestNodeSynchronizerFixture : public TestSynchronizerFixture {
public:
  static constexpr int max_int = std::numeric_limits<int>::max();

  void SetUp() override {
    TestSynchronizerFixture::SetUp();
    this->distribute();

    UInt nb_nodes = this->mesh->getNbNodes();

    node_data = std::make_unique<Array<int>>(nb_nodes);
    for (auto && data : enumerate(*node_data)) {
      auto n = std::get<0>(data);
      auto & d = std::get<1>(data);
      UInt gn = this->mesh->getNodeGlobalId(n);

      if (this->mesh->isMasterNode(n))
        d = gn;
      else if (this->mesh->isLocalNode(n))
        d = -gn;
      else if (this->mesh->isSlaveNode(n))
        d = max_int;
      else
        d = -max_int;
    }

    data_accessor = std::make_unique<DataAccessorTest>(*node_data);
  }

  void TearDown() override {
    data_accessor.reset(nullptr);
    node_data.reset(nullptr);
  }

  void checkData() {
    for (auto && data : enumerate(*this->node_data)) {
      auto n = std::get<0>(data);
      auto & d = std::get<1>(data);
      UInt gn = this->mesh->getNodeGlobalId(n);

      if (this->mesh->isMasterNode(n))
        EXPECT_EQ(d, gn);

      else if (this->mesh->isLocalNode(n))
        EXPECT_EQ(d, -gn);

      else if (this->mesh->isSlaveNode(n))
        EXPECT_EQ(d, gn);

      else

        EXPECT_EQ(d, -max_int);
    }
  }

protected:
  std::unique_ptr<Array<int>> node_data;
  std::unique_ptr<DataAccessorTest> data_accessor;
};

/* -------------------------------------------------------------------------- */
constexpr int TestNodeSynchronizerFixture::max_int;

/* -------------------------------------------------------------------------- */
TEST_F(TestNodeSynchronizerFixture, SynchroneOnce) {
  auto & synchronizer = this->mesh->getNodeSynchronizer();
  synchronizer.synchronizeOnce(*this->data_accessor, SynchronizationTag::_test);
  this->checkData();
}

/* -------------------------------------------------------------------------- */
TEST_F(TestNodeSynchronizerFixture, Synchrone) {
  auto & node_synchronizer = this->mesh->getNodeSynchronizer();
  node_synchronizer.synchronize(*this->data_accessor,
                                SynchronizationTag::_test);
  this->checkData();
}

/* -------------------------------------------------------------------------- */
TEST_F(TestNodeSynchronizerFixture, Asynchrone) {
  auto & synchronizer = this->mesh->getNodeSynchronizer();
  synchronizer.asynchronousSynchronize(*this->data_accessor,
                                       SynchronizationTag::_test);

  std::random_device r;
  std::default_random_engine engine(r());
  std::uniform_int_distribution<int> uniform_dist(10, 100);
  std::chrono::microseconds delay(uniform_dist(engine));

  std::this_thread::sleep_for(delay);

  synchronizer.waitEndSynchronize(*this->data_accessor,
                                  SynchronizationTag::_test);
  this->checkData();
}

/* -------------------------------------------------------------------------- */
TEST_F(TestNodeSynchronizerFixture, Gather) {
  auto & synchronizer = this->mesh->getNodeSynchronizer();

  const auto & comm = akantu::Communicator::getStaticCommunicator();
  Int prank = comm.whoAmI();

  if (prank == 0) {
    Array<int> all_data(this->mesh->getNbGlobalNodes());
    synchronizer.gather(*(this->node_data), all_data);
    for (auto && data : enumerate(all_data)) {
      EXPECT_EQ(std::get<0>(data), std::abs(std::get<1>(data)));
    }
  } else {
    synchronizer.gather(*(this->node_data));
  }
}

/* -------------------------------------------------------------------------- */
TEST_F(TestNodeSynchronizerFixture, Scatter) {
  Array<int> local_data(this->mesh->getNbNodes(), 1, this->max_int);
  auto & synchronizer = this->mesh->getNodeSynchronizer();
  
  if (prank == 0) {
    Array<int> all_data(this->mesh->getNbGlobalNodes());
    for (auto && data : enumerate(all_data)) {
      std::get<1>(data) = std::get<0>(data);
    }
    synchronizer.scatter(local_data, all_data);
  } else {
    synchronizer.scatter(local_data);
  }

  for (auto && data : enumerate(local_data)) {
    auto && n = std::get<0>(data);
    auto && d = std::get<1>(data);
    UInt gn = this->mesh->getNodeGlobalId(n);
    if(this->mesh->isPureGhostNode(n)) {
      EXPECT_EQ(d, this->max_int);
    } else {
      EXPECT_EQ(d, gn);
    }
  }
}
