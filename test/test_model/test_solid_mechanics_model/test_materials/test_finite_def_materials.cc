/**
 * @file   test_finite_def_materials.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Fri Nov 17 2017
 * @date last modification: Tue Feb 20 2018
 *
 * @brief
 *
 * @section LICENSE
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
 */
/* -------------------------------------------------------------------------- */
#include "test_gtest_utils.hh"
#include "test_material_fixtures.hh"
/* -------------------------------------------------------------------------- */
#include <material_neohookean.hh>
#include <solid_mechanics_model.hh>
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <type_traits>
/* -------------------------------------------------------------------------- */

using namespace akantu;

using mat_types = ::testing::Types<Traits<MaterialNeohookean, 1>,
                                   Traits<MaterialNeohookean, 2>,
                                   Traits<MaterialNeohookean, 3>>;

/*****************************************************************/

template <> void FriendMaterial<MaterialNeohookean<3>>::testComputeStress() {
  AKANTU_TO_IMPLEMENT();
}

/*****************************************************************/
template <>
void FriendMaterial<MaterialNeohookean<3>>::testComputeTangentModuli() {
  AKANTU_TO_IMPLEMENT();
}

/*****************************************************************/

template <> void FriendMaterial<MaterialNeohookean<3>>::testEnergyDensity() {
  AKANTU_TO_IMPLEMENT();
}

/*****************************************************************/

namespace {

template <typename T>
class TestFiniteDefMaterialFixture : public ::TestMaterialFixture<T> {};

TYPED_TEST_SUITE(TestFiniteDefMaterialFixture, mat_types);

TYPED_TEST(TestFiniteDefMaterialFixture, DISABLED_ComputeStress) {
  this->material->testComputeStress();
}
TYPED_TEST(TestFiniteDefMaterialFixture, DISABLED_EnergyDensity) {
  this->material->testEnergyDensity();
}
TYPED_TEST(TestFiniteDefMaterialFixture, DISABLED_ComputeTangentModuli) {
  this->material->testComputeTangentModuli();
}
TYPED_TEST(TestFiniteDefMaterialFixture, DISABLED_DefComputeCelerity) {
  this->material->testCelerity();
}
} // namespace
/*****************************************************************/
