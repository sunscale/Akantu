/**
 * @file   test_solid_mechanics_model_kinetic_energy.cc
 *
 * @author Tobias Brink <tobias.brink@epfl.ch>
 *
 * @date creation: Fri Nov 17 2017
 * @date last modification: Fri Jan 26 2018
 *
 * @brief  test kinetic energy
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
 * This test uses a linear elastic material with density = 1, Young's modulus =
 * 1, and Poisson's ratio = 0 and imposes a uniform velocity of 1. The volume of
 * the mesh is 1 and thus we have a mass of 1 and therefore a kinetic energy of
 * 0.5*m*v² = 0.5. The kind of constitutive law should not matter for this test,
 * so we use linear elastic. We perform 5 timesteps and check the solution every
 * time.
 *
 */
/* -------------------------------------------------------------------------- */
#include "../test_solid_mechanics_model_fixture.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

namespace {

void test_body(SolidMechanicsModel & model, AnalysisMethod analysis_method) {
  const auto spatial_dimension = model.getSpatialDimension();

  getStaticParser().parse("test_solid_mechanics_model_"
                          "kinetic_energy_material.dat");

  model.initFull(_analysis_method = analysis_method);
  model.assembleMassLumped();

  /// impose initial velocity of 1, it should remain constant
  auto & velo = model.getVelocity();

  for (auto && velov : make_view(velo, spatial_dimension)) {
    velov(_x) = 1;
  }

  /// set up timestep
  auto time_step = model.getStableTimeStep() * 0.8;
  model.setTimeStep(time_step);

  /// run five times and look at the kinetic energy
  for (uint i = 0; i < 5; ++i) {

    /// make a step
    model.solveStep();

    /// compare energy to analytical solution
    const Real E_ref = 0.5;
    auto E_kin = model.getEnergy("kinetic");

    EXPECT_NEAR(E_ref, E_kin, 1e-8);
  }
}

TYPED_TEST(TestSMMFixture, KineticEnergyImplicit) {
  test_body(*(this->model), _implicit_dynamic);
}

TYPED_TEST(TestSMMFixture, KineticEnergyExplicit) {
  test_body(*(this->model), _explicit_lumped_mass);
}

} // namespace
