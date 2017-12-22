/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#if defined(AKANTU_USE_PYBIND11)
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
    //py::initialize_interpreter();
#endif
  }
  // Override this to define how to tear down the environment.
  void TearDown() override {
    ::akantu::finalize();
#if defined(AKANTU_USE_PYBIND11)
    //py::finalize_interpreter();
#endif
  }

protected:
  int & argc;
  char **& argv;
};
}

int main(int argc, char ** argv) {
#if defined(AKANTU_USE_PYBIND11)
  py::scoped_interpreter guard{};
#endif

  ::testing::InitGoogleTest(&argc, argv);
  ::testing::AddGlobalTestEnvironment(new AkaEnvironment(argc, argv));

  return RUN_ALL_TESTS();
}
