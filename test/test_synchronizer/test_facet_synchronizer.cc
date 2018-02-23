/**
 * @file   test_facet_synchronizer.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 *
 * @brief  Facet synchronizer test
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */
/* -------------------------------------------------------------------------- */
#include "test_synchronizers_fixture.hh"
#include "test_data_accessor.hh"
/* -------------------------------------------------------------------------- */
#include "element_synchronizer.hh"
/* -------------------------------------------------------------------------- */
#include <random>
#include <chrono>
#include <thread>
/* -------------------------------------------------------------------------- */

class TestFacetSynchronizerFixture : public TestSynchronizerFixture {
public:
  void SetUp() override {
    TestSynchronizerFixture::SetUp();
    this->distribute();

    this->mesh->initMeshFacets();

    /// compute barycenter for each element
    barycenters = std::make_unique<ElementTypeMapArray<Real>>("barycenters", "", 0);
    this->initBarycenters(*barycenters, this->mesh->getMeshFacets());

    test_accessor = std::make_unique<TestAccessor>(*this->mesh, *this->barycenters);
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
  synchronizer.synchronizeOnce(*this->test_accessor, _gst_test);
}

/* -------------------------------------------------------------------------- */
TEST_F(TestFacetSynchronizerFixture, Synchrone) {
  auto & synchronizer = this->mesh->getMeshFacets().getElementSynchronizer();
  synchronizer.synchronize(*this->test_accessor, _gst_test);
}

/* -------------------------------------------------------------------------- */
TEST_F(TestFacetSynchronizerFixture, Asynchrone) {
  auto & synchronizer = this->mesh->getMeshFacets().getElementSynchronizer();
  synchronizer.asynchronousSynchronize(*this->test_accessor, _gst_test);

  std::random_device r;
  std::default_random_engine engine(r());
  std::uniform_int_distribution<int> uniform_dist(10, 100);
  std::chrono::microseconds delay(uniform_dist(engine));

  std::this_thread::sleep_for(delay);

  synchronizer.waitEndSynchronize(*this->test_accessor, _gst_test);
}
