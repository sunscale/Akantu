#include "aka_common.hh"
#include <gtest/gtest.h>

namespace {
class AkaEnvironment : public ::testing::Environment {
public:
  AkaEnvironment(int & argc, char **& argv) : argc(argc), argv(argv) {}
  // Override this to define how to set up the environment.
  void SetUp() override { ::akantu::initialize(argc, argv); }
  // Override this to define how to tear down the environment.
  void TearDown() override {  }

protected:
  int & argc;
  char **& argv;
};
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::AddGlobalTestEnvironment(new AkaEnvironment(argc, argv));

  return RUN_ALL_TESTS();
}
