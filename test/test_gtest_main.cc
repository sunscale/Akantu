/**
 * @file   test_gtest_main.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Nov 09 2017
 * @date last modification: Tue Jan 16 2018
 *
 * @brief  Main function for gtest tests
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
#include "communicator.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#if defined(AKANTU_TEST_USE_PYBIND11)
#include <pybind11/embed.h>
namespace py = pybind11;
#endif
/* -------------------------------------------------------------------------- */

namespace {
class AkaEnvironment : public ::testing::Environment {
public:
  AkaEnvironment(int & argc, char **& argv) : argc(argc), argv(argv) {}
  // Override this to define how to set up the environment.
  void SetUp() override {
    ::akantu::initialize(argc, argv);

#if defined(AKANTU_USE_PYBIND11)
// py::initialize_interpreter();
#endif
  }
  // Override this to define how to tear down the environment.
  void TearDown() override {
    ::akantu::finalize();
#if defined(AKANTU_USE_PYBIND11)
// py::finalize_interpreter();
#endif
  }

protected:
  int & argc;
  char **& argv;
};
} // namespace

int main(int argc, char ** argv) {
#if defined(AKANTU_TEST_USE_PYBIND11)
  py::scoped_interpreter guard{};
#endif

  ::testing::InitGoogleTest(&argc, argv);
  ::testing::AddGlobalTestEnvironment(new AkaEnvironment(argc, argv));

  ::testing::TestEventListeners & listeners =
      ::testing::UnitTest::GetInstance()->listeners();
  if (::akantu::Communicator::getStaticCommunicator().whoAmI() != 0) {
    delete listeners.Release(listeners.default_result_printer());
  }

  return RUN_ALL_TESTS();
}
