/* -------------------------------------------------------------------------- */
#include "patch_test_linear_heat_transfer_fixture.hh"
/* -------------------------------------------------------------------------- */

TYPED_TEST(TestPatchTestHTMLinear, Explicit) {
  this->initModel(_implicit_dynamic, "heat_transfer_input.dat");

  const auto & coordinates = this->mesh->getNodes();
  auto & temperature = this->model->getTemperature();
  // set the position of all nodes to the static solution
  for (auto && tuple : zip(make_view(coordinates, this->dim),
                           make_view(temperature, 1))) {
    this->setLinearDOF(std::get<1>(tuple), std::get<0>(tuple));
  }

  for (UInt s = 0; s < 100; ++s) {
    this->model->solveStep();
  }

  this->checkAll();
}
