/**
 * @file   test_material_thermal.cc
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Wed Aug 04 2010
 * @date last modification: Mon Jan 29 2018
 *
 * @brief  Test the thermal material
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
 */

/* -------------------------------------------------------------------------- */
#include "material_thermal.hh"
#include "solid_mechanics_model.hh"
#include "test_material_fixtures.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <type_traits>
/* -------------------------------------------------------------------------- */

using namespace akantu;

using mat_types =
    ::testing::Types<Traits<MaterialThermal, 1>, Traits<MaterialThermal, 2>,
                     Traits<MaterialThermal, 3>>;

/* -------------------------------------------------------------------------- */
template <> void FriendMaterial<MaterialThermal<3>>::testComputeStress() {
  Real E = 1.;
  Real nu = .3;
  Real alpha = 2;
  setParam("E", E);
  setParam("nu", nu);
  setParam("alpha", alpha);

  Real deltaT = 1;
  Real sigma = 0;
  this->computeStressOnQuad(sigma, deltaT);
  Real solution = -E / (1 - 2 * nu) * alpha * deltaT;
  auto error = std::abs(sigma - solution);
  ASSERT_NEAR(error, 0, 1e-14);
}

template <> void FriendMaterial<MaterialThermal<2>>::testComputeStress() {
  Real E = 1.;
  Real nu = .3;
  Real alpha = 2;
  setParam("E", E);
  setParam("nu", nu);
  setParam("alpha", alpha);

  Real deltaT = 1;
  Real sigma = 0;
  this->computeStressOnQuad(sigma, deltaT);
  Real solution = -E / (1 - 2 * nu) * alpha * deltaT;
  auto error = std::abs(sigma - solution);
  ASSERT_NEAR(error, 0, 1e-14);
}

template <> void FriendMaterial<MaterialThermal<1>>::testComputeStress() {
  Real E = 1.;
  Real nu = .3;
  Real alpha = 2;
  setParam("E", E);
  setParam("nu", nu);
  setParam("alpha", alpha);

  Real deltaT = 1;
  Real sigma = 0;
  this->computeStressOnQuad(sigma, deltaT);
  Real solution = -E * alpha * deltaT;
  auto error = std::abs(sigma - solution);
  ASSERT_NEAR(error, 0, 1e-14);
}

namespace {

template <typename T>
class TestMaterialThermalFixture : public ::TestMaterialFixture<T> {};

TYPED_TEST_SUITE(TestMaterialThermalFixture, mat_types, );

TYPED_TEST(TestMaterialThermalFixture, ThermalComputeStress) {
  this->material->testComputeStress();
}
} // namespace
