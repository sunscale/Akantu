/**
 * @file   test_solid_mechanics_model_work_dynamics.cc
 *
 * @author Tobias Brink <tobias.brink@epfl.ch>
 *
 * @date creation: Tue Dec 15 2017
 * @date last modification: Dec Nov 15 2017
 *
 * @brief  test work in dynamic simulations
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
 * @section LICENSE
 *
 * Copyright (©)  2017 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "test_solid_mechanics_model_fixture.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

namespace {

void test_body(SolidMechanicsModel & model, Mesh & mesh,
               AnalysisMethod analysis_method, UInt steps_needed) {
  const auto spatial_dimension = model.getSpatialDimension();

  getStaticParser().parse("test_solid_mechanics_model_"
                          "work_material.dat");

  /// model initialization
  model.initFull(_analysis_method = analysis_method);
  model.assembleMassLumped();

  /// Create a node group for Neumann BCs.
  auto & apply_force_grp = mesh.createNodeGroup("apply_force");
  auto & fixed_grp = mesh.createNodeGroup("fixed");

  const auto & pos = mesh.getNodes();
  auto & lower = mesh.getLowerBounds();
  auto & upper = mesh.getUpperBounds();

  UInt i = 0;
  for (auto && posv : make_view(pos, spatial_dimension)) {
    if (posv(_x) > upper(_x) - 1e-6) {
      apply_force_grp.add(i);
    } else if (posv(_x) < lower(_x) + 1e-6)  {
      fixed_grp.add(i);
    }
    ++i;
  }

  mesh.createElementGroupFromNodeGroup("el_apply_force",
                                       "apply_force",
                                       spatial_dimension - 1);
  mesh.createElementGroupFromNodeGroup("el_fixed",
                                       "fixed",
                                       spatial_dimension - 1);

  /// set up timestep
  auto time_step = model.getStableTimeStep() * 0.1;
  model.setTimeStep(time_step);

  /// Do the sim
  std::vector<Real> displacements{0.00, 0.01, -0.01};
  for (auto && u : displacements) {
    model.applyBC(BC::Dirichlet::FixedValue(u, _x), "el_fixed");

    Vector<Real> surface_traction(spatial_dimension);
    surface_traction(_x) = 0.5;
    if (spatial_dimension == 1) { //TODO: this is a hack to work
                                  //      around non-implemented
                                  //      BC::Neumann::FromTraction for 1D
      auto & force = model.getForce();
      for (auto && pair : zip(make_view(pos, spatial_dimension),
                              make_view(force, spatial_dimension))) {
        auto & posv = std::get<0>(pair);
        auto & forcev = std::get<1>(pair);
        if (posv(_x) > upper(_x) - 1e-6) {
          forcev(_x) = surface_traction(_x);
        }
      }
    } else {
      model.applyBC(BC::Neumann::FromTraction(surface_traction),
                    "el_apply_force");
    }

    // First, "equilibrate" a bit to get a reference state of total
    // energy and work. This is needed when we have a Dirichlet with
    // finite displacement on one side.
    for (UInt i = 0; i < 25; ++i) {
      model.solveStep();
    }
    // Again, work reported by Akantu is infinitesimal (dW) and we
    // need to integrate a while to get a decent value.
    double Etot0 = model.getEnergy("potential") + model.getEnergy("kinetic");
    double W = 0.0;
    for (UInt i = 0; i < steps_needed; ++i) {
      /// Solve.
      model.solveStep();

      const auto dW = model.getEnergy("external work");

      W += dW;
    }
    // Finally check.
    const auto Epot = model.getEnergy("potential");
    const auto Ekin = model.getEnergy("kinetic");

    EXPECT_NEAR(W, Ekin + Epot - Etot0, 5e-2); // Sadly not very exact
                                               // for such a coarse mesh.
  }
}

/* TODO: this is currently disabled for terrible results and performance    
TYPED_TEST(TestSMMFixtureBar, WorkImplicit) {
  test_body(*(this->model), *(this->mesh), _implicit_dynamic, 500);
}
*/

TYPED_TEST(TestSMMFixtureBar, WorkExplicit) {
  test_body(*(this->model), *(this->mesh), _explicit_lumped_mass, 200);
}

} // namespace
