/**
 * @file   test_voigt_helper.cc
 *
 * @author Nicolas Richart
 *
 * @date creation  mer nov 13 2019
 *
 * @brief unit tests for VoigtHelper
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
#include "test_gtest_utils.hh"
/* -------------------------------------------------------------------------- */
#include <aka_voigthelper.hh>
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <unordered_map>
/* -------------------------------------------------------------------------- */

using namespace akantu;

template <class Dim_v> class VoigtHelperFixture : public ::testing::Test {
protected:
  using voigt_h = VoigtHelper<Dim_v::value>;
  constexpr static UInt dim = Dim_v::value;

  VoigtHelperFixture() {
    switch (this->dim) {
    case 1: {
      indices.push_back({0, 0});
      matrix = Matrix<Real>{{10}};
      vector = Vector<Real>{10};
      vector_factor = Vector<Real>{10};
      break;
    }
    case 2: {
      indices.push_back({0, 0});
      indices.push_back({1, 1});
      indices.push_back({0, 1});

      matrix = Matrix<Real>{{10, 33}, {0, 56}};
      vector = Vector<Real>{10, 56, 33};
      vector_factor = Vector<Real>{10, 56, 2 * 33};
      break;
    }
    case 3: {
      indices.push_back({0, 0});
      indices.push_back({1, 1});
      indices.push_back({2, 2});
      indices.push_back({1, 2});
      indices.push_back({0, 2});
      indices.push_back({0, 1});

      matrix = Matrix<Real>{{10, 33, 20}, {0, 56, 27}, {0, 0, 98}};
      vector = Vector<Real>{10, 56, 98, 27, 20, 33};
      vector_factor = Vector<Real>{10, 56, 98, 2 * 27, 2 * 20, 2 * 33};
      break;
    }
    }
  }

  void SetUp() override {}

  std::vector<std::pair<UInt, UInt>> indices;
  Matrix<Real> matrix;
  Vector<Real> vector;
  Vector<Real> vector_factor;
};

template <UInt dim>
using spatial_dimension_t = std::integral_constant<UInt, dim>;

using TestTypes =
    ::testing::Types<spatial_dimension_t<1>, spatial_dimension_t<2>,
                     spatial_dimension_t<3>>;
TYPED_TEST_SUITE(VoigtHelperFixture, TestTypes);

TYPED_TEST(VoigtHelperFixture, Size) {
  using voigt_h = typename TestFixture::voigt_h;
  switch (this->dim) {
  case 1:
    EXPECT_EQ(voigt_h::size, 1);
    break;
  case 2:
    EXPECT_EQ(voigt_h::size, 3);
    break;
  case 3:
    EXPECT_EQ(voigt_h::size, 6);
    break;
  }
}

TYPED_TEST(VoigtHelperFixture, Indicies) {
  using voigt_h = typename TestFixture::voigt_h;

  for (UInt I = 0; I < voigt_h::size; ++I) {
    EXPECT_EQ(this->indices[I].first, voigt_h::vec[I][0]);
    EXPECT_EQ(this->indices[I].second, voigt_h::vec[I][1]);
  }
}

TYPED_TEST(VoigtHelperFixture, Factors) {
  using voigt_h = typename TestFixture::voigt_h;
  for (UInt I = 0; I < voigt_h::size; ++I) {
    if (I < this->dim) {
      EXPECT_EQ(voigt_h::factors[I], 1);
    } else {
      EXPECT_EQ(voigt_h::factors[I], 2);
    }
  }
}

TYPED_TEST(VoigtHelperFixture, MatrixToVoight) {
  using voigt_h = typename TestFixture::voigt_h;

  auto voigt = voigt_h::matrixToVoigt(this->matrix);

  for (UInt I = 0; I < voigt_h::size; ++I) {
    EXPECT_EQ(voigt(I), this->vector(I));
  }
}

TYPED_TEST(VoigtHelperFixture, MatrixToVoightFactors) {
  using voigt_h = typename TestFixture::voigt_h;

  auto voigt = voigt_h::matrixToVoigtWithFactors(this->matrix);

  for (UInt I = 0; I < voigt_h::size; ++I) {
    EXPECT_EQ(voigt(I), this->vector_factor(I));
  }
}

TYPED_TEST(VoigtHelperFixture, VoightToMatrix) {
  using voigt_h = typename TestFixture::voigt_h;

  auto matrix = voigt_h::voigtToMatrix(this->vector);

  for (UInt i = 0; i < this->dim; ++i) {
    for (UInt j = 0; j < this->dim; ++j) {
      EXPECT_EQ(matrix(i, j), this->matrix(std::min(i, j), std::max(i, j)));
    }
  }
}
