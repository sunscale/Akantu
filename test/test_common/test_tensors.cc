/**
 * @file   test_tensors.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Nov 14 2017
 * @date last modification: Mon Jan 22 2018
 *
 * @brief  test the tensors types
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
 */

/* -------------------------------------------------------------------------- */
#include "aka_array.hh"
#include "aka_iterators.hh"
#include "aka_types.hh"
/* -------------------------------------------------------------------------- */
#include <cstdlib>
#include <gtest/gtest.h>
#include <memory>
/* -------------------------------------------------------------------------- */

using namespace akantu;

namespace {

/* -------------------------------------------------------------------------- */
class TensorConstructorFixture : public ::testing::Test {
public:
  void SetUp() override {
    for (auto & r : reference) {
      r = rand(); // google-test seeds srand()
    }
  }
  void TearDown() override {}

  template <typename V> void compareToRef(const V & v) {
    for (int i = 0; i < size_; ++i) {
      EXPECT_DOUBLE_EQ(reference[i], v.storage()[i]);
    }
  }

protected:
  const int size_{24};
  const std::array<int, 2> mat_size{{4, 6}};
  // const std::array<int, 3> tens3_size{{4, 2, 3}};
  std::array<double, 24> reference;
};

/* -------------------------------------------------------------------------- */
class TensorFixture : public TensorConstructorFixture {
public:
  TensorFixture()
      : vref(reference.data(), size_),
        mref(reference.data(), mat_size[0], mat_size[1]) {}

protected:
  Vector<double> vref;
  Matrix<double> mref;
};

/* -------------------------------------------------------------------------- */
// Vector ----------------------------------------------------------------------
TEST_F(TensorConstructorFixture, VectorDefaultConstruct) {
  Vector<double> v;
  EXPECT_EQ(0, v.size());
  EXPECT_EQ(nullptr, v.storage());
  EXPECT_EQ(false, v.isWrapped());
}

TEST_F(TensorConstructorFixture, VectorConstruct1) {
  double r = rand();
  Vector<double> v(size_, r);
  EXPECT_EQ(size_, v.size());
  EXPECT_EQ(false, v.isWrapped());

  for (int i = 0; i < size_; ++i) {
    EXPECT_DOUBLE_EQ(r, v(i));
    EXPECT_DOUBLE_EQ(r, v[i]);
  }
}

TEST_F(TensorConstructorFixture, VectorConstructWrapped) {
  Vector<double> v(reference.data(), size_);
  EXPECT_EQ(size_, v.size());
  EXPECT_EQ(true, v.isWrapped());

  for (int i = 0; i < size_; ++i) {
    EXPECT_DOUBLE_EQ(reference[i], v(i));
    EXPECT_DOUBLE_EQ(reference[i], v[i]);
  }
}

TEST_F(TensorConstructorFixture, VectorConstructInitializer) {
  Vector<double> v{0., 1., 2., 3., 4., 5.};
  EXPECT_EQ(6, v.size());
  EXPECT_EQ(false, v.isWrapped());

  for (int i = 0; i < 6; ++i) {
    EXPECT_DOUBLE_EQ(i, v(i));
  }
}

TEST_F(TensorConstructorFixture, VectorConstructCopy1) {
  Vector<double> vref(reference.data(), reference.size());
  Vector<double> v(vref);
  EXPECT_EQ(size_, v.size());
  EXPECT_EQ(false, v.isWrapped());
  compareToRef(v);
}

TEST_F(TensorConstructorFixture, VectorConstructCopy2) {
  Vector<double> vref(reference.data(), reference.size());
  Vector<double> v(vref, false);
  EXPECT_EQ(size_, v.size());
  EXPECT_EQ(true, v.isWrapped());
  compareToRef(v);
}

TEST_F(TensorConstructorFixture, VectorConstructProxy1) {
  VectorProxy<double> vref(reference.data(), reference.size());
  EXPECT_EQ(size_, vref.size());
  compareToRef(vref);

  Vector<double> v(vref);
  EXPECT_EQ(size_, v.size());
  EXPECT_EQ(true, v.isWrapped());
  compareToRef(v);
}

TEST_F(TensorConstructorFixture, VectorConstructProxy2) {
  Vector<double> vref(reference.data(), reference.size());
  VectorProxy<double> v(vref);
  EXPECT_EQ(size_, v.size());
  compareToRef(v);
}

/* -------------------------------------------------------------------------- */
TEST_F(TensorFixture, VectorEqual) {
  Vector<double> v;

  v = vref;
  compareToRef(v);

  EXPECT_EQ(size_, v.size());
  EXPECT_EQ(false, v.isWrapped());
}

TEST_F(TensorFixture, VectorEqualProxy) {
  VectorProxy<double> vref_proxy(vref);
  Vector<double> v;

  v = vref;
  compareToRef(v);

  EXPECT_EQ(size_, v.size());
  EXPECT_EQ(false, v.isWrapped());
}

TEST_F(TensorFixture, VectorEqualProxy2) {
  Vector<double> v_store(size_, 0.);
  VectorProxy<double> v(v_store);

  v = vref;
  compareToRef(v);
  compareToRef(v_store);
}

/* -------------------------------------------------------------------------- */
TEST_F(TensorFixture, VectorSet) {
  Vector<double> v(vref);
  compareToRef(v);

  double r = rand();
  v.set(r);

  for (int i = 0; i < size_; ++i)
    EXPECT_DOUBLE_EQ(r, v[i]);
}

TEST_F(TensorFixture, VectorClear) {
  Vector<double> v(vref);
  compareToRef(v);

  v.zero();

  for (int i = 0; i < size_; ++i)
    EXPECT_DOUBLE_EQ(0, v[i]);
}

/* -------------------------------------------------------------------------- */
TEST_F(TensorFixture, VectorDivide) {
  Vector<double> v;
  double r = rand();
  v = vref / r;

  for (int i = 0; i < size_; ++i)
    EXPECT_DOUBLE_EQ(reference[i] / r, v[i]);
}

TEST_F(TensorFixture, VectorMultiply1) {
  Vector<double> v;
  double r = rand();
  v = vref * r;

  for (int i = 0; i < size_; ++i)
    EXPECT_DOUBLE_EQ(reference[i] * r, v[i]);
}

TEST_F(TensorFixture, VectorMultiply2) {
  Vector<double> v;
  double r = rand();
  v = r * vref;

  for (int i = 0; i < size_; ++i)
    EXPECT_DOUBLE_EQ(reference[i] * r, v[i]);
}

TEST_F(TensorFixture, VectorAddition) {
  Vector<double> v;
  v = vref + vref;

  for (int i = 0; i < size_; ++i)
    EXPECT_DOUBLE_EQ(reference[i] * 2., v[i]);
}

TEST_F(TensorFixture, VectorSubstract) {
  Vector<double> v;
  v = vref - vref;

  for (int i = 0; i < size_; ++i)
    EXPECT_DOUBLE_EQ(0., v[i]);
}

TEST_F(TensorFixture, VectorDivideEqual) {
  Vector<double> v(vref);
  double r = rand();
  v /= r;

  for (int i = 0; i < size_; ++i)
    EXPECT_DOUBLE_EQ(reference[i] / r, v[i]);
}

TEST_F(TensorFixture, VectorMultiplyEqual1) {
  Vector<double> v(vref);
  double r = rand();
  v *= r;

  for (int i = 0; i < size_; ++i)
    EXPECT_DOUBLE_EQ(reference[i] * r, v[i]);
}

TEST_F(TensorFixture, VectorMultiplyEqual2) {
  Vector<double> v(vref);
  v *= v;

  for (int i = 0; i < size_; ++i)
    EXPECT_DOUBLE_EQ(reference[i] * reference[i], v[i]);
}

TEST_F(TensorFixture, VectorAdditionEqual) {
  Vector<double> v(vref);
  v += vref;

  for (int i = 0; i < size_; ++i)
    EXPECT_DOUBLE_EQ(reference[i] * 2., v[i]);
}

TEST_F(TensorFixture, VectorSubstractEqual) {
  Vector<double> v(vref);
  v -= vref;

  for (int i = 0; i < size_; ++i)
    EXPECT_DOUBLE_EQ(0., v[i]);
}

/* -------------------------------------------------------------------------- */
// Matrix ----------------------------------------------------------------------

TEST_F(TensorConstructorFixture, MatrixDefaultConstruct) {
  Matrix<double> m;
  EXPECT_EQ(0, m.size());
  EXPECT_EQ(0, m.rows());
  EXPECT_EQ(0, m.cols());
  EXPECT_EQ(nullptr, m.storage());
  EXPECT_EQ(false, m.isWrapped());
}

TEST_F(TensorConstructorFixture, MatrixConstruct1) {
  double r = rand();
  Matrix<double> m(mat_size[0], mat_size[1], r);
  EXPECT_EQ(size_, m.size());
  EXPECT_EQ(mat_size[0], m.rows());
  EXPECT_EQ(mat_size[1], m.cols());
  EXPECT_EQ(false, m.isWrapped());

  for (int i = 0; i < mat_size[0]; ++i) {
    for (int j = 0; j < mat_size[1]; ++j) {
      EXPECT_EQ(r, m(i, j));
      EXPECT_EQ(r, m[i + j * mat_size[0]]);
    }
  }
}

TEST_F(TensorConstructorFixture, MatrixConstructWrapped) {
  Matrix<double> m(reference.data(), mat_size[0], mat_size[1]);
  EXPECT_EQ(size_, m.size());
  EXPECT_EQ(mat_size[0], m.rows());
  EXPECT_EQ(mat_size[1], m.cols());
  EXPECT_EQ(true, m.isWrapped());

  for (int i = 0; i < mat_size[0]; ++i) {
    for (int j = 0; j < mat_size[1]; ++j) {
      EXPECT_DOUBLE_EQ(reference[i + j * mat_size[0]], m(i, j));
    }
  }
  compareToRef(m);
}

TEST_F(TensorConstructorFixture, MatrixConstructInitializer) {
  Matrix<double> m{{0., 1., 2.}, {3., 4., 5.}};
  EXPECT_EQ(6, m.size());
  EXPECT_EQ(2, m.rows());
  EXPECT_EQ(3, m.cols());

  EXPECT_EQ(false, m.isWrapped());

  int c = 0;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 3; ++j, ++c) {
      EXPECT_DOUBLE_EQ(c, m(i, j));
    }
  }
}

TEST_F(TensorConstructorFixture, MatrixConstructCopy1) {
  Matrix<double> mref(reference.data(), mat_size[0], mat_size[1]);
  Matrix<double> m(mref);
  EXPECT_EQ(size_, m.size());
  EXPECT_EQ(mat_size[0], m.rows());
  EXPECT_EQ(mat_size[1], m.cols());
  EXPECT_EQ(false, m.isWrapped());
  compareToRef(m);
}

TEST_F(TensorConstructorFixture, MatrixConstructCopy2) {
  Matrix<double> mref(reference.data(), mat_size[0], mat_size[1]);
  Matrix<double> m(mref);
  EXPECT_EQ(size_, m.size());
  EXPECT_EQ(mat_size[0], m.rows());
  EXPECT_EQ(mat_size[1], m.cols());
  EXPECT_EQ(false, m.isWrapped());
  compareToRef(m);
}

TEST_F(TensorConstructorFixture, MatrixConstructProxy1) {
  MatrixProxy<double> mref(reference.data(), mat_size[0], mat_size[1]);
  EXPECT_EQ(size_, mref.size());
  EXPECT_EQ(mat_size[0], mref.size(0));
  EXPECT_EQ(mat_size[1], mref.size(1));
  compareToRef(mref);

  Matrix<double> m(mref);
  EXPECT_EQ(size_, m.size());
  EXPECT_EQ(mat_size[0], m.rows());
  EXPECT_EQ(mat_size[1], m.cols());
  EXPECT_EQ(true, m.isWrapped());
  compareToRef(m);
}

TEST_F(TensorConstructorFixture, MatrixConstructProxy2) {
  Matrix<double> mref(reference.data(), mat_size[0], mat_size[1]);
  MatrixProxy<double> m(mref);
  EXPECT_EQ(size_, m.size());
  EXPECT_EQ(mat_size[0], m.size(0));
  EXPECT_EQ(mat_size[1], m.size(1));
  compareToRef(m);
}

/* -------------------------------------------------------------------------- */
TEST_F(TensorFixture, MatrixEqual) {
  Matrix<double> m;

  m = mref;
  compareToRef(m);

  EXPECT_EQ(size_, m.size());
  EXPECT_EQ(mat_size[0], m.rows());
  EXPECT_EQ(mat_size[1], m.cols());
  EXPECT_EQ(false, m.isWrapped());
}

TEST_F(TensorFixture, MatrixEqualProxy1) {
  MatrixProxy<double> mref_proxy(mref);
  Matrix<double> m;

  m = mref;
  compareToRef(m);

  EXPECT_EQ(size_, m.size());
  EXPECT_EQ(mat_size[0], m.rows());
  EXPECT_EQ(mat_size[1], m.cols());
  EXPECT_EQ(false, m.isWrapped());
}

TEST_F(TensorFixture, MatrixEqualProxy2) {
  Matrix<double> m_store(mat_size[0], mat_size[1], 0.);
  MatrixProxy<double> m(m_store);

  m = mref;
  compareToRef(m);
  compareToRef(m_store);
}

TEST_F(TensorFixture, MatrixEqualSlice) {
  Matrix<double> m(mat_size[0], mat_size[1], 0.);

  for (unsigned int i = 0; i < m.cols(); ++i)
    m(i) = Vector<Real>(mref(i));

  compareToRef(m);
}

/* -------------------------------------------------------------------------- */
TEST_F(TensorFixture, MatrixSet) {
  Matrix<double> m(mref);
  compareToRef(m);

  double r = rand();
  m.set(r);

  for (int i = 0; i < size_; ++i)
    EXPECT_DOUBLE_EQ(r, m[i]);
}

TEST_F(TensorFixture, MatrixClear) {
  Matrix<double> m(mref);
  compareToRef(m);

  m.zero();

  for (int i = 0; i < size_; ++i)
    EXPECT_DOUBLE_EQ(0, m[i]);
}

/* -------------------------------------------------------------------------- */
TEST_F(TensorFixture, MatrixDivide) {
  Matrix<double> m;
  double r = rand();
  m = mref / r;

  for (int i = 0; i < size_; ++i)
    EXPECT_DOUBLE_EQ(reference[i] / r, m[i]);
}

TEST_F(TensorFixture, MatrixMultiply1) {
  Matrix<double> m;
  double r = rand();
  m = mref * r;

  for (int i = 0; i < size_; ++i)
    EXPECT_DOUBLE_EQ(reference[i] * r, m[i]);
}

TEST_F(TensorFixture, MatrixMultiply2) {
  Matrix<double> m;
  double r = rand();
  m = r * mref;

  for (int i = 0; i < size_; ++i)
    EXPECT_DOUBLE_EQ(reference[i] * r, m[i]);
}

TEST_F(TensorFixture, MatrixAddition) {
  Matrix<double> m;
  m = mref + mref;

  for (int i = 0; i < size_; ++i)
    EXPECT_DOUBLE_EQ(reference[i] * 2., m[i]);
}

TEST_F(TensorFixture, MatrixSubstract) {
  Matrix<double> m;
  m = mref - mref;

  for (int i = 0; i < size_; ++i)
    EXPECT_DOUBLE_EQ(0., m[i]);
}

TEST_F(TensorFixture, MatrixDivideEqual) {
  Matrix<double> m(mref);
  double r = rand();
  m /= r;

  for (int i = 0; i < size_; ++i)
    EXPECT_DOUBLE_EQ(reference[i] / r, m[i]);
}

TEST_F(TensorFixture, MatrixMultiplyEqual1) {
  Matrix<double> m(mref);
  double r = rand();
  m *= r;

  for (int i = 0; i < size_; ++i)
    EXPECT_DOUBLE_EQ(reference[i] * r, m[i]);
}

TEST_F(TensorFixture, MatrixAdditionEqual) {
  Matrix<double> m(mref);
  m += mref;

  for (int i = 0; i < size_; ++i)
    EXPECT_DOUBLE_EQ(reference[i] * 2., m[i]);
}

TEST_F(TensorFixture, MatrixSubstractEqual) {
  Matrix<double> m(mref);
  m -= mref;

  for (int i = 0; i < size_; ++i)
    EXPECT_DOUBLE_EQ(0., m[i]);
}

TEST_F(TensorFixture, MatrixIterator) {
  Matrix<double> m(mref);

  UInt col_count = 0;
  for (auto && col : m) {
    Vector<Real> col_hand(m.storage() + col_count * m.rows(), m.rows());
    Vector<Real> col_wrap(col);

    auto comp = (col_wrap - col_hand).norm<L_inf>();
    EXPECT_DOUBLE_EQ(0., comp);
    ++col_count;
  }
}

TEST_F(TensorFixture, MatrixIteratorZip) {
  Matrix<double> m1(mref);
  Matrix<double> m2(mref);

  UInt col_count = 0;
  for (auto && col : zip(m1, m2)) {
    Vector<Real> col1(std::get<0>(col));
    Vector<Real> col2(std::get<1>(col));

    auto comp = (col1 - col2).norm<L_inf>();
    EXPECT_DOUBLE_EQ(0., comp);
    ++col_count;
  }
}

#if defined(AKANTU_USE_LAPACK)
TEST_F(TensorFixture, MatrixEigs) {
  Matrix<double> m{{0, 1, 0, 0}, {1., 0, 0, 0}, {0, 1, 0, 1}, {0, 0, 4, 0}};

  Matrix<double> eig_vects(4, 4);
  Vector<double> eigs(4);
  m.eig(eigs, eig_vects);

  Vector<double> eigs_ref{2, 1., -1., -2};
  auto lambda_v = m * eig_vects;

  for (int i = 0; i < 4; ++i) {
    EXPECT_NEAR(eigs_ref(i), eigs(i), 1e-14);
    for (int j = 0; j < 4; ++j) {
      EXPECT_NEAR(lambda_v(i)(j), eigs(i) * eig_vects(i)(j), 1e-14);
    }
  }
}
#endif

/* -------------------------------------------------------------------------- */

} // namespace
