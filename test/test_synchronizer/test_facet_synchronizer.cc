/**
 * @file   test_facet_synchronizer.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Nov 05 2014
 * @date last modification: Fri Jan 26 2018
 *
 * @brief  Facet synchronizer test
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "test_data_accessor.hh"
#include "test_synchronizers_fixture.hh"
/* -------------------------------------------------------------------------- */
#include "element_synchronizer.hh"
/* -------------------------------------------------------------------------- */
#include <chrono>
#include <random>
#include <thread>
/* -------------------------------------------------------------------------- */

class TestFacetSynchronizerFixture : public TestSynchronizerFixture {
public:
  void SetUp() override {
    TestSynchronizerFixture::SetUp();
    this->distribute();

    this->mesh->initMeshFacets();

    /// compute barycenter for each element
    barycenters =
        std::make_unique<ElementTypeMapArray<Real>>("barycenters", "", 0);
    this->initBarycenters(*barycenters, this->mesh->getMeshFacets());

    test_accessor =
        std::make_unique<TestAccessor>(*this->mesh, *this->barycenters);
  }

  void TearDown() override {
    barycenters.reset(nullptr);
    test_accessor.reset(nullptr);
  }

protected:
  std::unique_ptr<ElementTypeMapArray<Real>> barycenters;
  std::unique_ptr<TestAccessor> test_accessor;
};

/* -------------------------------------------------------------------------- */
TEST_F(TestFacetSynchronizerFixture, SynchroneOnce) {
  auto & synchronizer = this->mesh->getMeshFacets().getElementSynchronizer();
  synchronizer.synchronizeOnce(*this->test_accessor, SynchronizationTag::_test);
}

/* -------------------------------------------------------------------------- */
TEST_F(TestFacetSynchronizerFixture, Synchrone) {
  auto & synchronizer = this->mesh->getMeshFacets().getElementSynchronizer();
  synchronizer.synchronize(*this->test_accessor, SynchronizationTag::_test);
}

/* -------------------------------------------------------------------------- */
TEST_F(TestFacetSynchronizerFixture, Asynchrone) {
  auto & synchronizer = this->mesh->getMeshFacets().getElementSynchronizer();
  synchronizer.asynchronousSynchronize(*this->test_accessor,
                                       SynchronizationTag::_test);

  std::random_device r;
  std::default_random_engine engine(r());
  std::uniform_int_distribution<int> uniform_dist(10, 100);
  std::chrono::microseconds delay(uniform_dist(engine));

  std::this_thread::sleep_for(delay);

  synchronizer.waitEndSynchronize(*this->test_accessor,
                                  SynchronizationTag::_test);
}
