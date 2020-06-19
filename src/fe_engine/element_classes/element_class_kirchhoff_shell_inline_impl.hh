/**
 * @file   element_class_kirchhoff_shell_inline_impl.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Damien Spielmann <damien.spielmann@epfl.ch>
 *
 * @date creation: Fri Jul 04 2014
 * @date last modification: Wed Feb 21 2018
 *
 * @brief  Element class Kirchhoff Shell
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "element_class_structural.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_ELEMENT_CLASS_KIRCHHOFF_SHELL_INLINE_IMPL_HH__
#define __AKANTU_ELEMENT_CLASS_KIRCHHOFF_SHELL_INLINE_IMPL_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
AKANTU_DEFINE_STRUCTURAL_INTERPOLATION_TYPE_PROPERTY(
    _itp_discrete_kirchhoff_triangle_18, _itp_lagrange_triangle_3, 6, 6, 21);
AKANTU_DEFINE_STRUCTURAL_ELEMENT_CLASS_PROPERTY(
    _discrete_kirchhoff_triangle_18, _gt_triangle_3,
    _itp_discrete_kirchhoff_triangle_18, _triangle_3, _ek_structural, 3,
    _git_triangle, 2);

/* -------------------------------------------------------------------------- */
namespace detail {
  inline void computeBasisChangeMatrix(Matrix<Real> & P,
                                       const Matrix<Real> & X) {
    Vector<Real> X1 = X(0);
    Vector<Real> X2 = X(1);
    Vector<Real> X3 = X(2);

    Vector<Real> a1 = X2 - X1;
    Vector<Real> a2 = X3 - X1;

    a1.normalize();
    Vector<Real> e3 = a1.crossProduct(a2);
    e3.normalize();
    Vector<Real> e2 = e3.crossProduct(a1);

    P(0) = a1;
    P(1) = e2;
    P(2) = e3;
    P = P.transpose();
  }
} // namespace detail

/* -------------------------------------------------------------------------- */
template <>
inline void
ElementClass<_discrete_kirchhoff_triangle_18>::computeRotationMatrix(
    Matrix<Real> & R, const Matrix<Real> & X, const Vector<Real> &) {
  auto dim = X.rows();
  Matrix<Real> P(dim, dim);
  detail::computeBasisChangeMatrix(P, X);

  R.clear();
  for (UInt i = 0; i < dim; ++i)
    for (UInt j = 0; j < dim; ++j)
      R(i + dim, j + dim) = R(i, j) = P(i, j);
}

/* -------------------------------------------------------------------------- */
template <>
inline void
InterpolationElement<_itp_discrete_kirchhoff_triangle_18>::computeShapes(
    const Vector<Real> & /*natural_coords*/,
    const Matrix<Real> & /*real_coord*/, Matrix<Real> & /*N*/) {}

/* -------------------------------------------------------------------------- */
template <>
inline void
InterpolationElement<_itp_discrete_kirchhoff_triangle_18>::computeDNDS(
    const Vector<Real> & natural_coords, const Matrix<Real> & real_coordinates,
    Matrix<Real> & B) {

  auto dim = real_coordinates.cols();
  Matrix<Real> P(dim, dim);
  detail::computeBasisChangeMatrix(P, real_coordinates);
  auto X = P * real_coordinates;
  Vector<Real> X1 = X(0);
  Vector<Real> X2 = X(1);
  Vector<Real> X3 = X(2);

  std::array<Vector<Real>, 3> A = {X2 - X1, X3 - X2, X1 - X3};
  std::array<Real, 3> L, C, S;

  // Setting all last coordinates to 0
  std::for_each(A.begin(), A.end(), [](auto & a) { a(2) = 0; });
  // Computing lengths
  std::transform(A.begin(), A.end(), L.begin(),
                 [](auto & a) { return a.template norm<L_2>(); });
  // Computing cosines
  std::transform(A.begin(), A.end(), L.begin(), C.begin(),
                 [](auto & a, auto & l) { return a(0) / l; });
  // Computing sines
  std::transform(A.begin(), A.end(), L.begin(), S.begin(),
                 [](auto & a, auto & l) { return a(1) / l; });

  // Natural coordinates
  Real xi = natural_coords(0);
  Real eta = natural_coords(1);

  // Derivative of quadratic interpolation functions
  Matrix<Real> dP = {{4 * (1 - 2 * xi - eta), 4 * eta, -4 * eta},
                     {-4 * xi, 4 * xi, 4 * (1 - xi - 2 * eta)}};

  Matrix<Real> dNx1 = {
      {3. / 2 * (dP(0, 0) * C[0] / L[0] - dP(0, 2) * C[2] / L[2]),
       3. / 2 * (dP(0, 1) * C[1] / L[1] - dP(0, 0) * C[0] / L[0]),
       3. / 2 * (dP(0, 2) * C[2] / L[2] - dP(0, 1) * C[1] / L[1])},
      {3. / 2 * (dP(1, 0) * C[0] / L[0] - dP(1, 2) * C[2] / L[2]),
       3. / 2 * (dP(1, 1) * C[1] / L[1] - dP(1, 0) * C[0] / L[0]),
       3. / 2 * (dP(1, 2) * C[2] / L[2] - dP(1, 1) * C[1] / L[1])}};
  Matrix<Real> dNx2 = {
      // clang-format off
      {-1 - 3. / 4 * (dP(0, 0) * C[0] * C[0] + dP(0, 2) * C[2] * C[2]),
	1 - 3. / 4 * (dP(0, 1) * C[1] * C[1] + dP(0, 0) * C[0] * C[0]),
	  - 3. / 4 * (dP(0, 2) * C[2] * C[2] + dP(0, 1) * C[1] * C[1])},
      {-1 - 3. / 4 * (dP(1, 0) * C[0] * C[0] + dP(1, 2) * C[2] * C[2]),
	  - 3. / 4 * (dP(1, 1) * C[1] * C[1] + dP(1, 0) * C[0] * C[0]),
	1 - 3. / 4 * (dP(1, 2) * C[2] * C[2] + dP(1, 1) * C[1] * C[1])}};
  // clang-format on
  Matrix<Real> dNx3 = {
      {-3. / 4 * (dP(0, 0) * C[0] * S[0] + dP(0, 2) * C[2] * S[2]),
       -3. / 4 * (dP(0, 1) * C[1] * S[1] + dP(0, 0) * C[0] * S[0]),
       -3. / 4 * (dP(0, 2) * C[2] * S[2] + dP(0, 1) * C[1] * S[1])},
      {-3. / 4 * (dP(1, 0) * C[0] * S[0] + dP(1, 2) * C[2] * S[2]),
       -3. / 4 * (dP(1, 1) * C[1] * S[1] + dP(1, 0) * C[0] * S[0]),
       -3. / 4 * (dP(1, 2) * C[2] * S[2] + dP(1, 1) * C[1] * S[1])}};
  Matrix<Real> dNy1 = {
      {3. / 2 * (dP(0, 0) * S[0] / L[0] - dP(0, 2) * S[2] / L[2]),
       3. / 2 * (dP(0, 1) * S[1] / L[1] - dP(0, 0) * S[0] / L[0]),
       3. / 2 * (dP(0, 2) * S[2] / L[2] - dP(0, 1) * S[1] / L[1])},
      {3. / 2 * (dP(1, 0) * S[0] / L[0] - dP(1, 2) * S[2] / L[2]),
       3. / 2 * (dP(1, 1) * S[1] / L[1] - dP(1, 0) * S[0] / L[0]),
       3. / 2 * (dP(1, 2) * S[2] / L[2] - dP(1, 1) * S[1] / L[1])}};
  Matrix<Real> dNy2 = dNx3;
  Matrix<Real> dNy3 = {
      // clang-format off
      {-1 - 3. / 4 * (dP(0, 0) * S[0] * S[0] + dP(0, 2) * S[2] * S[2]),
	1 - 3. / 4 * (dP(0, 1) * S[1] * S[1] + dP(0, 0) * S[0] * S[0]),
	  - 3. / 4 * (dP(0, 2) * S[2] * S[2] + dP(0, 1) * S[1] * S[1])},
      {-1 - 3. / 4 * (dP(1, 0) * S[0] * S[0] + dP(1, 2) * S[2] * S[2]),
	  - 3. / 4 * (dP(1, 1) * S[1] * S[1] + dP(1, 0) * S[0] * S[0]),
	1 - 3. / 4 * (dP(1, 2) * S[2] * S[2] + dP(1, 1) * S[1] * S[1])}};
  // clang-format on

  // Derivative of linear (membrane mode) functions
  Matrix<Real> dNm(2, 3);
  InterpolationElement<_itp_lagrange_triangle_3, _itk_lagrangian>::computeDNDS(
      natural_coords, dNm);

  UInt i = 0;
  for (const Matrix<Real> & mat : {dNm, dNx1, dNx2, dNx3, dNy1, dNy2, dNy3}) {
    B.block(mat, 0, i);
    i += mat.cols();
  }
}
/* -------------------------------------------------------------------------- */
template <>
inline void
InterpolationElement<_itp_discrete_kirchhoff_triangle_18,
                     _itk_structural>::arrangeInVoigt(const Matrix<Real> & dnds,
                                                      Matrix<Real> & B) {
  Matrix<Real> dNm(2, 3), dNx1(2, 3), dNx2(2, 3), dNx3(2, 3), dNy1(2, 3),
      dNy2(2, 3), dNy3(2, 3);
  UInt i = 0;
  for (Matrix<Real> * mat : {&dNm, &dNx1, &dNx2, &dNx3, &dNy1, &dNy2, &dNy3}) {
    *mat = dnds.block(0, i, 2, 3);
    i += mat->cols();
  }

  for (UInt i = 0; i < 3; ++i) {
    // clang-format off
    Matrix<Real> Bm = {{dNm(0, i), 0,         0, 0, 0, 0},
                       {0,         dNm(1, i), 0, 0, 0, 0},
                       {dNm(1, i), dNm(0, i), 0, 0, 0, 0}};
    Matrix<Real> Bf = {{0, 0, dNx1(0, i),              -dNx3(0, i),              dNx2(0, i),              0},
                       {0, 0, dNy1(1, i),              -dNy3(1, i),              dNy2(1, i),              0},
                       {0, 0, dNx1(1, i) + dNy1(0, i), -dNx3(1, i) - dNy3(0, i), dNx2(1, i) + dNy2(0, i), 0}};
    // clang-format on
    B.block(Bm, 0, i * 6);
    B.block(Bf, 3, i * 6);
  }
}

} // namespace akantu

#endif /* __AKANTU_ELEMENT_CLASS_KIRCHHOFF_SHELL_INLINE_IMPL_HH__ */
