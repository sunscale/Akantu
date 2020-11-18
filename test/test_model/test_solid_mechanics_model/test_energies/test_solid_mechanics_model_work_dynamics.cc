/**
 * @file   test_solid_mechanics_model_work_dynamics.cc
 *
 * @author Tobias Brink <tobias.brink@epfl.ch>
 *
 * @date creation: Fri Dec 15 2017
 * @date last modification: Fri Jan 26 2018
 *
 * @brief  test work in dynamic simulations
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
 * @section description
 *
 * Assuming that the kinetic energy and the potential energy of the
 * linear elastic material are bug free, the work in a dynamic
 * simulation must equal the change in internal energy (first law of
 * thermodynamics). Work in dynamics is an infinitesimal work Fds,
 * thus we need to integrate it and compare at the end. In this test,
 * we use one Dirichlet boundary condition (with u = 0.0, 0.01, and
 * -0.01) and one Neumann boundary condition for F on the opposite
 * side. Then we do a few steps to get reference energies for work and
 * internal energy. After more steps, the change in both work and
 * internal energy must be equal.
 *
 */

/* -------------------------------------------------------------------------- */
#include "../test_solid_mechanics_model_fixture.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

namespace {

template <typename type_>
class TestSMMFixtureWorkDynamic : public TestSMMFixture<type_> {
public:
  void SetUp() override {
    this->mesh_file =
        "../../patch_tests/data/bar" + std::to_string(this->type) + ".msh";
    TestSMMFixture<type_>::SetUp();

    getStaticParser().parse("test_solid_mechanics_model_"
                            "work_material.dat");

    /// model initialization
    this->model->initFull();

    /// Create a node group for Neumann BCs.
    auto & apply_force_grp = this->mesh->createNodeGroup("apply_force");
    auto & fixed_grp = this->mesh->createNodeGroup("fixed");

    const auto & pos = this->mesh->getNodes();
    const auto & lower = this->mesh->getLowerBounds();
    const auto & upper = this->mesh->getUpperBounds();

    UInt i = 0;
    for (auto && posv : make_view(pos, this->spatial_dimension)) {
      if (posv(_x) > upper(_x) - 1e-6) {
        apply_force_grp.add(i);
      } else if (posv(_x) < lower(_x) + 1e-6) {
        fixed_grp.add(i);
      }
      ++i;
    }

    this->mesh->createElementGroupFromNodeGroup("el_apply_force", "apply_force",
                                                this->spatial_dimension - 1);
    this->mesh->createElementGroupFromNodeGroup("el_fixed", "fixed",
                                                this->spatial_dimension - 1);

    Vector<Real> surface_traction(this->spatial_dimension);
    surface_traction(_x) = 0.5;
    if (this->spatial_dimension == 1) {
      // TODO: this is a hack to work
      //      around non-implemented
      //      BC::Neumann::FromTraction for 1D
      auto & force = this->model->getExternalForce();

      for (auto && pair : zip(make_view(pos, this->spatial_dimension),
                              make_view(force, this->spatial_dimension))) {
        auto & posv = std::get<0>(pair);
        auto & forcev = std::get<1>(pair);
        if (posv(_x) > upper(_x) - 1e-6) {
          forcev(_x) = surface_traction(_x);
        }
      }
    } else {
      this->model->applyBC(BC::Neumann::FromTraction(surface_traction),
                           "el_apply_force");
    }

    /// set up timestep
    auto time_step = this->model->getStableTimeStep() * 0.1;
    this->model->setTimeStep(time_step);
  }
};

TYPED_TEST_SUITE(TestSMMFixtureWorkDynamic, gtest_element_types, );

/* TODO: this is currently disabled for terrible results and performance
TYPED_TEST(TestSMMFixtureBar, WorkImplicit) {
  test_body(*(this->model), *(this->mesh), _implicit_dynamic, 500);
}
*/
// model.assembleMassLumped();
TYPED_TEST(TestSMMFixtureWorkDynamic, WorkExplicit) {
  /// Do the sim
  std::vector<Real> displacements{0.00, 0.01, -0.01};

  for (auto && u : displacements) {
    this->model->applyBC(BC::Dirichlet::FixedValue(u, _x), "el_fixed");

    // First, "equilibrate" a bit to get a reference state of total
    // energy and work. This is needed when we have a Dirichlet with
    // finite displacement on one side.
    for (UInt i = 0; i < 25; ++i) {
      this->model->solveStep();
    }
    // Again, work reported by Akantu is infinitesimal (dW) and we
    // need to integrate a while to get a decent value.
    double Etot0 =
        this->model->getEnergy("potential") + this->model->getEnergy("kinetic");
    double W = 0.0;
    for (UInt i = 0; i < 200; ++i) {
      /// Solve.
      this->model->solveStep();

      const auto dW = this->model->getEnergy("external work");

      W += dW;
    }
    // Finally check.
    const auto Epot = this->model->getEnergy("potential");
    const auto Ekin = this->model->getEnergy("kinetic");

    EXPECT_NEAR(W, Ekin + Epot - Etot0, 5e-2);
    // Sadly not very exact for such a coarse mesh.
  }
}
} // namespace
