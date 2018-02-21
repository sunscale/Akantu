/**
 * @file   test_material_cohesive_linear.cc
 *
 * @author Nicolas Richart
 *
 * @date creation  Tue Feb 20 2018
 *
 * @brief Test material cohesive linear
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
#include "test_material_cohesive_fixture.hh"
/* -------------------------------------------------------------------------- */
#include "material_cohesive_linear.hh"
/* -------------------------------------------------------------------------- */

template <UInt dim>
struct TestMaterialCohesiveLinear
    : public TestMaterialCohesive<MaterialCohesiveLinear, dim> {

  TestMaterialCohesiveLinear(SolidMechanicsModel & model)
      : TestMaterialCohesive<MaterialCohesiveLinear, dim>(model) {}

  void SetUp() override {
    this->is_extrinsic = true;

    this->beta = 1.;
    this->G_c = 10.;
    this->sigma_c_ = 1e6;
    this->delta_c_ = 2. * this->G_c / this->sigma_c_;
  }

  void resetInternal() override {
    normal_opening = Vector<Real>(dim, 0.);
    tangential_opening = Vector<Real>(dim, 0.);
    contact_traction = Vector<Real>(dim, 0.);
    contact_opening = Vector<Real>(dim, 0.);
  }

  void computeTractions(const Array<Real> & openings,
                        const Vector<Real> & normal,
                        Array<Real> & tractions) override {
    for (auto && data :
         zip(make_view(openings, dim), make_view(tractions, dim))) {
      const auto & opening = std::get<0>(data);
      auto & traction = std::get<1>(data);

      this->computeTractionOnQuad(
          traction, opening, normal, delta_max, this->delta_c_,
          this->insertion_stress_, this->sigma_c_, normal_opening,
          tangential_opening, normal_opening_norm, tangential_opening_norm,
          damage, penetration, contact_traction, contact_opening);
    }
  }

public:
  Real delta_c_{0};

  Real delta_max{0.};
  Real normal_opening_norm{0};
  Real tangential_opening_norm{0};
  Real damage{0};
  bool penetration{false};

  Vector<Real> normal_opening;
  Vector<Real> tangential_opening;

  Vector<Real> contact_traction;
  Vector<Real> contact_opening;
};

template <typename dim_>
using TestMaterialCohesiveLinearFixture =
    TestMaterialCohesiveFixture<TestMaterialCohesiveLinear, dim_>;

using types = gtest_list_t<TestAllDimensions>;

TYPED_TEST_CASE(TestMaterialCohesiveLinearFixture, types);

TYPED_TEST(TestMaterialCohesiveLinearFixture, ModeI) {
  auto normal = this->applyModIOpening(this->material->delta_c_, 100);

  auto delta_c = this->material->delta_c_;
  auto sigma_c = this->material->sigma_c_;

  for (auto && data : zip(make_view(*this->openings, this->dim),
                          make_view(*this->tractions, this->dim))) {
    const auto & opening = std::get<0>(data);
    auto & traction = std::get<1>(data);

    auto delta = opening.dot(normal);
    auto T = traction.dot(normal);

    auto T_expected = sigma_c * (delta_c - delta) / delta_c;

    EXPECT_NEAR(T_expected, T, 1e-9);
  }
}

// TYPED_TEST(TestMaterialCohesiveLinearFixture, ModeII) {
//   auto normal = this->applyModIIOpening(this->material->delta_c_, 100);

//   auto delta_c = this->material->delta_c_;
//   auto sigma_c = this->material->sigma_c_;

//   for (auto && data : zip(make_view(*this->openings, this->dim),
//                           make_view(*this->tractions, this->dim))) {
//     const auto & opening = std::get<0>(data);
//     auto & traction = std::get<1>(data);

//     auto delta = opening.dot(normal);

//     auto T = ;

//     auto T_expected = sigma_c * (delta_c - delta) / delta_c;

//     EXPECT_NEAR(T_expected, T, 1e-9);
//   }
// }
