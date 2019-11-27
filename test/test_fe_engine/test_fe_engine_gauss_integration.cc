/**
 * @file   test_fe_engine_gauss_integration.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue May 24 2016
 * @date last modification: Mon Feb 19 2018
 *
 * @brief  test integration on elements, this test consider that mesh is a cube
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
#include "test_fe_engine_fixture.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <iostream>
/* -------------------------------------------------------------------------- */
using namespace akantu;

namespace {
/* -------------------------------------------------------------------------- */
template <size_t t> using degree_t = std::integral_constant<size_t, t>;

/* -------------------------------------------------------------------------- */
using TestDegreeTypes = std::tuple<degree_t<0>, degree_t<1>, degree_t<2>,
                                   degree_t<3>, degree_t<4>, degree_t<5>>;

std::array<Polynomial<5>, 3> global_polys{
    {{0.40062394, 0.13703225, 0.51731446, 0.87830084, 0.5410543, 0.71842292},
     {0.41861835, 0.11080576, 0.49874043, 0.49077504, 0.85073835, 0.66259755},
     {0.92620845, 0.7503478, 0.62962232, 0.31662719, 0.64069644, 0.30878135}}};

template <typename T>
class TestGaussIntegrationFixture
    : public TestFEMFixture<std::tuple_element_t<0, T>> {
protected:
  using parent = TestFEMFixture<std::tuple_element_t<0, T>>;
  static constexpr size_t degree{std::tuple_element_t<1, T>::value};

public:
  TestGaussIntegrationFixture() : integration_points_pos(0, parent::dim) {}

  void SetUp() override {
    parent::SetUp();
    this->fem->initShapeFunctions();

    auto integration_points =
             this->fem->getIntegrator().template getIntegrationPoints <
             parent::type,
         degree == 0 ? 1 : degree > ();

    nb_integration_points = integration_points.cols();

    auto shapes_size = ElementClass<parent::type>::getShapeSize();
    Array<Real> shapes(0, shapes_size);
    this->fem->getShapeFunctions()
        .template computeShapesOnIntegrationPoints<parent::type>(
            this->mesh->getNodes(), integration_points, shapes, _not_ghost);

    auto vect_size = this->nb_integration_points * this->nb_element;
    integration_points_pos.resize(vect_size);
    this->fem->getShapeFunctions()
        .template interpolateOnIntegrationPoints<parent::type>(
            this->mesh->getNodes(), integration_points_pos, this->dim, shapes);

    for (size_t d = 0; d < this->dim; ++d) {
      polys[d] = global_polys[d].extract(degree);
    }
  }

  void testIntegrate() {
    std::stringstream sstr;
    sstr << this->type << ":" << this->degree;
    SCOPED_TRACE(sstr.str().c_str());

    auto vect_size = this->nb_integration_points * this->nb_element;
    Array<Real> polynomial(vect_size);
    size_t dim = parent::dim;

    for (size_t d = 0; d < dim; ++d) {
      auto poly = this->polys[d];
      for (auto && pair :
           zip(polynomial, make_view(this->integration_points_pos, dim))) {
        auto && p = std::get<0>(pair);
        auto & x = std::get<1>(pair);
        p = poly(x(d));
      }

      auto res =
          this->fem->getIntegrator()
              .template integrate<parent::type, (degree == 0 ? 1 : degree)>(
                  polynomial);
      auto expect = poly.integrate(this->lower(d), this->upper(d));

      for (size_t o = 0; o < dim; ++o) {
        if (o == d)
          continue;
        expect *= this->upper(d) - this->lower(d);
      }

      EXPECT_NEAR(expect, res, 5e-14);
    }
  }

protected:
  UInt nb_integration_points;
  std::array<Array<Real>, parent::dim> polynomial;
  Array<Real> integration_points_pos;
  std::array<Polynomial<5>, 3> polys;
};

template <typename T> constexpr size_t TestGaussIntegrationFixture<T>::degree;

/* -------------------------------------------------------------------------- */
/* Tests                                                                      */
/* -------------------------------------------------------------------------- */
TYPED_TEST_SUITE_P(TestGaussIntegrationFixture);

TYPED_TEST_P(TestGaussIntegrationFixture, ArbitraryOrder) {
  this->testIntegrate();
}

REGISTER_TYPED_TEST_SUITE_P(TestGaussIntegrationFixture, ArbitraryOrder);

using TestTypes = gtest_list_t<
    tuple_split_t<50, cross_product_t<TestElementTypes, TestDegreeTypes>>>;

INSTANTIATE_TYPED_TEST_SUITE_P(Split1, TestGaussIntegrationFixture, TestTypes);

using TestTypesTail = gtest_list_t<
    tuple_split_tail_t<50, cross_product_t<TestElementTypes, TestDegreeTypes>>>;

INSTANTIATE_TYPED_TEST_SUITE_P(Split2, TestGaussIntegrationFixture,
                               TestTypesTail);
} // namespace
