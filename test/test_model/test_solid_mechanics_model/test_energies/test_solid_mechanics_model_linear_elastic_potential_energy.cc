/**
 * @file   test_solid_mechanics_model_linear_elastic_potential_energy.cc
 *
 * @author Tobias Brink <tobias.brink@epfl.ch>
 *
 * @date creation: Tue Nov 14 2017
 * @date last modification: Fri Jan 26 2018
 *
 * @brief  test potential energy of the linear elasticity model
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
 * This test uses a linear elastic material with density = 1, Young's
 * modulus = 1, and Poisson's ratio = 0 and applies a linear
 * displacement from 0 to ε in x direction. The resulting potential
 * energy density should be 0.5*Y*ε² = ε²/2. Since the mesh always has
 * a volume of 1, the energy density equals the total energy. We test
 * 3 different strains.
 *
 */
/* -------------------------------------------------------------------------- */
#include "../test_solid_mechanics_model_fixture.hh"
#include "sparse_matrix.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

namespace {

TYPED_TEST(TestSMMFixture, LinearElasticPotentialEnergy) {
  const auto spatial_dimension = this->spatial_dimension;

  this->initModel("test_solid_mechanics_model_linear_elastic_"
                  "potential_energy_material.dat",
                  _static);

  const auto & lower = this->mesh->getLowerBounds();
  const auto & upper = this->mesh->getUpperBounds();
  auto length = upper(_x) - lower(_x);

  const auto & pos = this->mesh->getNodes();
  auto & disp = this->model->getDisplacement();
  auto & boun = this->model->getBlockedDOFs();

  std::vector<Real> strains{0.0, 0.1, 0.2, 0.3};

  for (auto && eps : strains) {
    /// boundary conditions
    for (auto && pair : zip(make_view(pos, spatial_dimension),
                            make_view(disp, spatial_dimension),
                            make_view(boun, spatial_dimension))) {
      const auto & posv = std::get<0>(pair);
      auto & dispv = std::get<1>(pair);
      auto & bounv = std::get<2>(pair);
      auto reduced_x = (posv(_x) - lower(_x)) / length;
      dispv(_x) = reduced_x * eps;
      bounv(_x) = true;

      if (posv(_x) < (lower(_x) + 1e-6)) {
        if ((spatial_dimension > 1) and (posv(_y) < (lower(_y) + 1e-6))) {
          bounv(_y) = true;

          if ((spatial_dimension > 2) and (posv(_z) < (lower(_z) + 1e-6))) {
            bounv(_z) = true;
          }
        }
      }
    }

    if (this->dump_paraview) {
      this->model->dump();
    }
    /// "solve" a step (solution is imposed)
    try {
      this->model->solveStep();
    } catch (...) {
      const auto & A = this->model->getDOFManager().getMatrix("J");
      auto prank = this->mesh->getCommunicator().whoAmI();
      A.saveMatrix("solver_mumps" + std::to_string(prank) + ".mtx");
      throw;
    }

    if (this->dump_paraview) {
      const auto & A = this->model->getDOFManager().getMatrix("J");
      auto prank = this->mesh->getCommunicator().whoAmI();
      A.saveMatrix("solver_mumps" + std::to_string(prank) + ".mtx");
    }

    /// compare energy to analytical solution
    auto E_ref = 0.5 * eps * eps;
    auto E_pot = this->model->getEnergy("potential");

    EXPECT_NEAR(E_ref, E_pot, 1e-8);
  }
}

} // namespace
