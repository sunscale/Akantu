/* -------------------------------------------------------------------------- */
#include "mesh_utils.hh"
#include "non_linear_solver.hh"
#include "solid_mechanics_model.hh"
#include "test_gtest_utils.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <vector>
/* -------------------------------------------------------------------------- */

TYPED_TEST(TestPatchTestLinearElastic, Implicit) {
  this->model->initFull(_analysis_method = _static);

  this->applyBC();

  auto & solver = this->model->getNonLinearSolver();
  solver.set("max_iterations", 2);
  solver.set("threshold", 2e-4);
  solver.set("convergence_type", _scc_residual);

  this->model->solveStep();

  this->checkDisplacements();
  this->checkStrains();
  this->checkStresses();
}

TYPED_TEST(TestPatchTestLinearElastic, Explicit) {
  this->model->initFull(_analysis_method = _explicit_lumped_mass);

  this->applyBC();

  const auto & coordinates = this->mesh->getNodes();
  auto & displacement = this->model->getDisplacement();

  //set the position of all nodes to the static solution
  for (auto && tuple :
           zip(make_view(coordinates, this->dim), make_view(displacement, this->dim))) {
    this->setDisplacement(std::get<1>(tuple), std::get<0>(tuple));
  }

  for(UInt s = 0; s < 100; ++s) {
    this->model->solveStep();
  }

  auto ekin = this->model->getEnergy("kinetic");
  EXPECT_NEAR(0, ekin, 1e-16);

  this->checkDisplacements();
  this->checkStrains();
  this->checkStresses();
}
