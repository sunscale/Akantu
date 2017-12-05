/* -------------------------------------------------------------------------- */
#include "patch_test_linear_fixture.hh"
/* -------------------------------------------------------------------------- */

TYPED_TEST(TestPatchTestLinear, Implicit) {
  std::string filename = "material_check_stress_plane_stress.dat";
  if (this->plane_strain)
    filename = "material_check_stress_plane_strain.dat";

  this->initModel(_static, filename);

  auto & solver = this->model->getNonLinearSolver();
  solver.set("max_iterations", 2);
  solver.set("threshold", 2e-4);
  solver.set("convergence_type", _scc_residual);

  this->model->solveStep();

  this->checkDisplacements();
  this->checkStrains();
  this->checkStresses([&](const Matrix<Real> & pstrain) {
    Real nu = this->model->getMaterial(0).get("nu");
    Real E = this->model->getMaterial(0).get("E");

    auto strain = (pstrain + pstrain.transpose()) / 2.;
    auto trace = strain.trace();

    auto lambda = nu * E / ((1 + nu) * (1 - 2 * nu));
    auto mu = E / (2 * (1 + nu));

    if (not this->plane_strain) {
      lambda = nu * E / (1 - nu * nu);
    }

    decltype(strain) stress(this->dim, this->dim);

    if (this->dim == 1) {
      stress(0, 0) = E * strain(0, 0);
    } else {
      for (UInt i = 0; i < this->dim; ++i)
        for (UInt j = 0; j < this->dim; ++j)
          stress(i, j) = (i == j) * lambda * trace + 2 * mu * strain(i, j);
    }

    return stress;
  });
}
