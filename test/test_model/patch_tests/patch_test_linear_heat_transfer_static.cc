/* -------------------------------------------------------------------------- */
#include "patch_test_linear_heat_transfer_fixture.hh"
/* -------------------------------------------------------------------------- */
#include "non_linear_solver.hh"
/* -------------------------------------------------------------------------- */

TYPED_TEST(TestPatchTestHTMLinear, Static) {
  this->initModel(_static, "heat_transfer_input.dat");

  auto & solver = this->model->getNonLinearSolver();
  solver.set("max_iterations", 2);
  solver.set("threshold", 2e-4);
  solver.set("convergence_type", _scc_residual);

  this->model->solveStep();

  this->checkAll();
}
