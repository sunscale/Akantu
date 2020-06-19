/**
 * @file   aka_math_tmpl.hh
 *
 * @author Ramin Aghababaei <ramin.aghababaei@epfl.ch>
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Mathilde Radiguet <mathilde.radiguet@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
 * @author Peter Spijker <peter.spijker@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Aug 04 2010
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Implementation of the inline functions of the math toolkit
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
#include "aka_blas_lapack.hh"
#include "aka_math.hh"
/* -------------------------------------------------------------------------- */
#include <cmath>
#include <cstring>
#include <typeinfo>
/* -------------------------------------------------------------------------- */

namespace akantu {
/* -------------------------------------------------------------------------- */
inline void Math::matrix_vector(UInt im, UInt in, Real * A, Real * x, Real * y,
                                Real alpha) {
#ifdef AKANTU_USE_BLAS
  /// y = alpha*op(A)*x + beta*y
  char tran_A = 'N';
  int incx = 1;
  int incy = 1;
  double beta = 0.;
  int m = im;
  int n = in;

  aka_gemv(&tran_A, &m, &n, &alpha, A, &m, x, &incx, &beta, y, &incy);

#else
  memset(y, 0, im * sizeof(Real));
  for (UInt i = 0; i < im; ++i) {
    for (UInt j = 0; j < in; ++j) {
      y[i] += A[i + j * im] * x[j];
    }
    y[i] *= alpha;
  }
#endif
}

/* -------------------------------------------------------------------------- */
inline void Math::matrixt_vector(UInt im, UInt in, Real * A, Real * x, Real * y,
                                 Real alpha) {
#ifdef AKANTU_USE_BLAS
  /// y = alpha*op(A)*x + beta*y
  char tran_A = 'T';
  int incx = 1;
  int incy = 1;
  double beta = 0.;
  int m = im;
  int n = in;

  aka_gemv(&tran_A, &m, &n, &alpha, A, &m, x, &incx, &beta, y, &incy);
#else
  memset(y, 0, in * sizeof(Real));
  for (UInt i = 0; i < im; ++i) {
    for (UInt j = 0; j < in; ++j) {
      y[j] += A[j * im + i] * x[i];
    }
    y[i] *= alpha;
  }
#endif
}

/* -------------------------------------------------------------------------- */
inline void Math::matrix_matrix(UInt im, UInt in, UInt ik, Real * A, Real * B,
                                Real * C, Real alpha) {
#ifdef AKANTU_USE_BLAS
  ///  C := alpha*op(A)*op(B) + beta*C
  char trans_a = 'N';
  char trans_b = 'N';
  double beta = 0.;
  int m = im, n = in, k = ik;

  aka_gemm(&trans_a, &trans_b, &m, &n, &k, &alpha, A, &m, B, &k, &beta, C, &m);
#else
  memset(C, 0, im * in * sizeof(Real));
  for (UInt j = 0; j < in; ++j) {
    UInt _jb = j * ik;
    UInt _jc = j * im;
    for (UInt i = 0; i < im; ++i) {
      for (UInt l = 0; l < ik; ++l) {
        UInt _la = l * im;
        C[i + _jc] += A[i + _la] * B[l + _jb];
      }
      C[i + _jc] *= alpha;
    }
  }
#endif
}

/* -------------------------------------------------------------------------- */
inline void Math::matrixt_matrix(UInt im, UInt in, UInt ik, Real * A, Real * B,
                                 Real * C, Real alpha) {
#ifdef AKANTU_USE_BLAS
  ///  C := alpha*op(A)*op(B) + beta*C
  char trans_a = 'T';
  char trans_b = 'N';
  double beta = 0.;
  int m = im, n = in, k = ik;

  aka_gemm(&trans_a, &trans_b, &m, &n, &k, &alpha, A, &k, B, &k, &beta, C, &m);
#else
  memset(C, 0, im * in * sizeof(Real));
  for (UInt j = 0; j < in; ++j) {
    UInt _jc = j * im;
    UInt _jb = j * ik;
    for (UInt i = 0; i < im; ++i) {
      UInt _ia = i * ik;
      for (UInt l = 0; l < ik; ++l) {
        C[i + _jc] += A[l + _ia] * B[l + _jb];
      }
      C[i + _jc] *= alpha;
    }
  }
#endif
}

/* -------------------------------------------------------------------------- */
inline void Math::matrix_matrixt(UInt im, UInt in, UInt ik, Real * A, Real * B,
                                 Real * C, Real alpha) {
#ifdef AKANTU_USE_BLAS
  ///  C := alpha*op(A)*op(B) + beta*C
  char trans_a = 'N';
  char trans_b = 'T';
  double beta = 0.;
  int m = im, n = in, k = ik;

  aka_gemm(&trans_a, &trans_b, &m, &n, &k, &alpha, A, &m, B, &n, &beta, C, &m);
#else
  memset(C, 0, im * in * sizeof(Real));
  for (UInt j = 0; j < in; ++j) {
    UInt _jc = j * im;
    for (UInt i = 0; i < im; ++i) {
      for (UInt l = 0; l < ik; ++l) {
        UInt _la = l * im;
        UInt _lb = l * in;
        C[i + _jc] += A[i + _la] * B[j + _lb];
      }
      C[i + _jc] *= alpha;
    }
  }
#endif
}

/* -------------------------------------------------------------------------- */
inline void Math::matrixt_matrixt(UInt im, UInt in, UInt ik, Real * A, Real * B,
                                  Real * C, Real alpha) {
#ifdef AKANTU_USE_BLAS
  ///  C := alpha*op(A)*op(B) + beta*C
  char trans_a = 'T';
  char trans_b = 'T';
  double beta = 0.;
  int m = im, n = in, k = ik;

  aka_gemm(&trans_a, &trans_b, &m, &n, &k, &alpha, A, &k, B, &n, &beta, C, &m);
#else
  memset(C, 0, im * in * sizeof(Real));
  for (UInt j = 0; j < in; ++j) {
    UInt _jc = j * im;
    for (UInt i = 0; i < im; ++i) {
      UInt _ia = i * ik;
      for (UInt l = 0; l < ik; ++l) {
        UInt _lb = l * in;
        C[i + _jc] += A[l + _ia] * B[j + _lb];
      }
      C[i + _jc] *= alpha;
    }
  }
#endif
}

/* -------------------------------------------------------------------------- */
inline void Math::aXplusY(UInt n, Real alpha, Real * x, Real * y) {
#ifdef AKANTU_USE_BLAS
  ///  y := alpha x + y
  int incx = 1, incy = 1;
  aka_axpy(&n, &alpha, x, &incx, y, &incy);
#else
  for (UInt i = 0; i < n; ++i)
    *(y++) += alpha * *(x++);
#endif
}

/* -------------------------------------------------------------------------- */
inline Real Math::vectorDot(Real * v1, Real * v2, UInt in) {
#ifdef AKANTU_USE_BLAS
  ///  d := v1 . v2
  int incx = 1, incy = 1, n = in;
  Real d = aka_dot(&n, v1, &incx, v2, &incy);
#else
  Real d = 0;
  for (UInt i = 0; i < in; ++i) {
    d += v1[i] * v2[i];
  }
#endif
  return d;
}

/* -------------------------------------------------------------------------- */
template <bool tr_A, bool tr_B>
inline void Math::matMul(UInt m, UInt n, UInt k, Real alpha, Real * A, Real * B,
                         __attribute__((unused)) Real beta, Real * C) {
  if (tr_A) {
    if (tr_B)
      matrixt_matrixt(m, n, k, A, B, C, alpha);
    else
      matrixt_matrix(m, n, k, A, B, C, alpha);
  } else {
    if (tr_B)
      matrix_matrixt(m, n, k, A, B, C, alpha);
    else
      matrix_matrix(m, n, k, A, B, C, alpha);
  }
}

/* -------------------------------------------------------------------------- */
template <bool tr_A>
inline void Math::matVectMul(UInt m, UInt n, Real alpha, Real * A, Real * x,
                             __attribute__((unused)) Real beta, Real * y) {
  if (tr_A) {
    matrixt_vector(m, n, A, x, y, alpha);
  } else {
    matrix_vector(m, n, A, x, y, alpha);
  }
}

/* -------------------------------------------------------------------------- */
template <typename T> inline void Math::matrixEig(UInt n, T * A, T * d, T * V) {

  // Matrix  A is  row major,  so the  lapack function  in fortran  will process
  // A^t. Asking for the left eigenvectors of A^t will give the transposed right
  // eigenvectors of A so in the C++ code the right eigenvectors.
  char jobvr, jobvl;
  if (V != nullptr)
    jobvr = 'V'; // compute left  eigenvectors
  else
    jobvr = 'N'; // compute left  eigenvectors

  jobvl = 'N'; // compute right eigenvectors

  auto * di = new T[n]; // imaginary part of the eigenvalues

  int info;
  int N = n;

  T wkopt;
  int lwork = -1;
  // query and allocate the optimal workspace
  aka_geev<T>(&jobvl, &jobvr, &N, A, &N, d, di, nullptr, &N, V, &N, &wkopt,
              &lwork, &info);

  lwork = int(wkopt);
  auto * work = new T[lwork];
  // solve the eigenproblem
  aka_geev<T>(&jobvl, &jobvr, &N, A, &N, d, di, nullptr, &N, V, &N, work,
              &lwork, &info);

  AKANTU_DEBUG_ASSERT(
      info == 0,
      "Problem computing eigenvalues/vectors. DGEEV exited with the value "
          << info);

  delete[] work;
  delete[] di; // I hope for you that there was no complex eigenvalues !!!
}

/* -------------------------------------------------------------------------- */
inline void Math::matrix22_eigenvalues(Real * A, Real * Adiag) {
  /// d = determinant of Matrix A
  Real d = det2(A);
  /// b = trace of Matrix A
  Real b = A[0] + A[3];

  Real c = sqrt(b * b - 4 * d);
  Adiag[0] = .5 * (b + c);
  Adiag[1] = .5 * (b - c);
}

/* -------------------------------------------------------------------------- */
inline void Math::matrix33_eigenvalues(Real * A, Real * Adiag) {
  matrixEig(3, A, Adiag);
}

/* -------------------------------------------------------------------------- */
template <UInt dim> inline void Math::eigenvalues(Real * A, Real * d) {
  if (dim == 1) {
    d[0] = A[0];
  } else if (dim == 2) {
    matrix22_eigenvalues(A, d);
  }
  // else if(dim == 3) { matrix33_eigenvalues(A, d); }
  else
    matrixEig(dim, A, d);
}

/* -------------------------------------------------------------------------- */
inline Real Math::det2(const Real * mat) {
  return mat[0] * mat[3] - mat[1] * mat[2];
}

/* -------------------------------------------------------------------------- */
inline Real Math::det3(const Real * mat) {
  return mat[0] * (mat[4] * mat[8] - mat[7] * mat[5]) -
         mat[3] * (mat[1] * mat[8] - mat[7] * mat[2]) +
         mat[6] * (mat[1] * mat[5] - mat[4] * mat[2]);
}

/* -------------------------------------------------------------------------- */
template <UInt n> inline Real Math::det(const Real * mat) {
  if (n == 1)
    return *mat;
  else if (n == 2)
    return det2(mat);
  else if (n == 3)
    return det3(mat);
  else
    return det(n, mat);
}

/* -------------------------------------------------------------------------- */
template <typename T> inline T Math::det(UInt n, const T * A) {
  int N = n;
  int info;
  auto * ipiv = new int[N + 1];

  auto * LU = new T[N * N];
  std::copy(A, A + N * N, LU);

  // LU factorization of A
  aka_getrf(&N, &N, LU, &N, ipiv, &info);
  if (info > 0) {
    AKANTU_ERROR("Singular matrix - cannot factorize it (info: " << info
                                                                 << " )");
  }

  // det(A) = det(L) * det(U) = 1 * det(U) = product_i U_{ii}
  T det = 1.;
  for (int i = 0; i < N; ++i)
    det *= (2 * (ipiv[i] == i) - 1) * LU[i * n + i];

  delete[] ipiv;
  delete[] LU;
  return det;
}

/* -------------------------------------------------------------------------- */
inline void Math::normal2(const Real * vec, Real * normal) {
  normal[0] = vec[1];
  normal[1] = -vec[0];
  Math::normalize2(normal);
}

/* -------------------------------------------------------------------------- */
inline void Math::normal3(const Real * vec1, const Real * vec2, Real * normal) {
  Math::vectorProduct3(vec1, vec2, normal);
  Math::normalize3(normal);
}

/* -------------------------------------------------------------------------- */
inline void Math::normalize2(Real * vec) {
  Real norm = Math::norm2(vec);
  vec[0] /= norm;
  vec[1] /= norm;
}

/* -------------------------------------------------------------------------- */
inline void Math::normalize3(Real * vec) {
  Real norm = Math::norm3(vec);
  vec[0] /= norm;
  vec[1] /= norm;
  vec[2] /= norm;
}

/* -------------------------------------------------------------------------- */
inline Real Math::norm2(const Real * vec) {
  return sqrt(vec[0] * vec[0] + vec[1] * vec[1]);
}

/* -------------------------------------------------------------------------- */
inline Real Math::norm3(const Real * vec) {
  return sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}

/* -------------------------------------------------------------------------- */
inline Real Math::norm(UInt n, const Real * vec) {
  Real norm = 0.;
  for (UInt i = 0; i < n; ++i) {
    norm += vec[i] * vec[i];
  }
  return sqrt(norm);
}

/* -------------------------------------------------------------------------- */
inline void Math::inv2(const Real * mat, Real * inv) {
  Real det_mat = det2(mat);

  inv[0] = mat[3] / det_mat;
  inv[1] = -mat[1] / det_mat;
  inv[2] = -mat[2] / det_mat;
  inv[3] = mat[0] / det_mat;
}

/* -------------------------------------------------------------------------- */
inline void Math::inv3(const Real * mat, Real * inv) {
  Real det_mat = det3(mat);

  inv[0] = (mat[4] * mat[8] - mat[7] * mat[5]) / det_mat;
  inv[1] = (mat[2] * mat[7] - mat[8] * mat[1]) / det_mat;
  inv[2] = (mat[1] * mat[5] - mat[4] * mat[2]) / det_mat;
  inv[3] = (mat[5] * mat[6] - mat[8] * mat[3]) / det_mat;
  inv[4] = (mat[0] * mat[8] - mat[6] * mat[2]) / det_mat;
  inv[5] = (mat[2] * mat[3] - mat[5] * mat[0]) / det_mat;
  inv[6] = (mat[3] * mat[7] - mat[6] * mat[4]) / det_mat;
  inv[7] = (mat[1] * mat[6] - mat[7] * mat[0]) / det_mat;
  inv[8] = (mat[0] * mat[4] - mat[3] * mat[1]) / det_mat;
}

/* -------------------------------------------------------------------------- */
template <UInt n> inline void Math::inv(const Real * A, Real * Ainv) {
  if (n == 1)
    *Ainv = 1. / *A;
  else if (n == 2)
    inv2(A, Ainv);
  else if (n == 3)
    inv3(A, Ainv);
  else
    inv(n, A, Ainv);
}

/* -------------------------------------------------------------------------- */
template <typename T> inline void Math::inv(UInt n, const T * A, T * invA) {
  int N = n;
  int info;
  auto * ipiv = new int[N + 1];
  int lwork = N * N;
  auto * work = new T[lwork];

  std::copy(A, A + n * n, invA);

  aka_getrf(&N, &N, invA, &N, ipiv, &info);
  if (info > 0) {
    AKANTU_ERROR("Singular matrix - cannot factorize it (info: " << info
                                                                 << " )");
  }

  aka_getri(&N, invA, &N, ipiv, work, &lwork, &info);
  if (info != 0) {
    AKANTU_ERROR("Cannot invert the matrix (info: " << info << " )");
  }

  delete[] ipiv;
  delete[] work;
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline void Math::solve(UInt n, const T * A, T * x, const T * b) {
  int N = n;
  int info;
  auto * ipiv = new int[N];
  auto * lu_A = new T[N * N];

  std::copy(A, A + N * N, lu_A);

  aka_getrf(&N, &N, lu_A, &N, ipiv, &info);
  if (info > 0) {
    AKANTU_ERROR("Singular matrix - cannot factorize it (info: " << info
                                                                 << " )");
  }

  char trans = 'N';
  int nrhs = 1;

  std::copy(b, b + N, x);

  aka_getrs(&trans, &N, &nrhs, lu_A, &N, ipiv, x, &N, &info);
  if (info != 0) {
    AKANTU_ERROR("Cannot solve the system (info: " << info << " )");
  }

  delete[] ipiv;
  delete[] lu_A;
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
inline Real Math::matrixDoubleDot22(Real * A, Real * B) {
  Real d;
  d = A[0] * B[0] + A[1] * B[1] + A[2] * B[2] + A[3] * B[3];
  return d;
}

/* -------------------------------------------------------------------------- */
inline Real Math::matrixDoubleDot33(Real * A, Real * B) {
  Real d;
  d = A[0] * B[0] + A[1] * B[1] + A[2] * B[2] + A[3] * B[3] + A[4] * B[4] +
      A[5] * B[5] + A[6] * B[6] + A[7] * B[7] + A[8] * B[8];
  return d;
}

/* -------------------------------------------------------------------------- */
inline Real Math::matrixDoubleDot(UInt n, Real * A, Real * B) {
  Real d = 0.;
  for (UInt i = 0; i < n; ++i) {
    for (UInt j = 0; j < n; ++j) {
      d += A[i * n + j] * B[i * n + j];
    }
  }
  return d;
}

/* -------------------------------------------------------------------------- */
inline void Math::vectorProduct3(const Real * v1, const Real * v2, Real * res) {
  res[0] = v1[1] * v2[2] - v1[2] * v2[1];
  res[1] = v1[2] * v2[0] - v1[0] * v2[2];
  res[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

/* -------------------------------------------------------------------------- */
inline Real Math::vectorDot2(const Real * v1, const Real * v2) {
  return (v1[0] * v2[0] + v1[1] * v2[1]);
}

/* -------------------------------------------------------------------------- */
inline Real Math::vectorDot3(const Real * v1, const Real * v2) {
  return (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]);
}

/* -------------------------------------------------------------------------- */
inline Real Math::distance_2d(const Real * x, const Real * y) {
  return sqrt((y[0] - x[0]) * (y[0] - x[0]) + (y[1] - x[1]) * (y[1] - x[1]));
}

/* -------------------------------------------------------------------------- */
inline Real Math::triangle_inradius(const Real * coord1, const Real * coord2,
                                    const Real * coord3) {
  /**
   * @f{eqnarray*}{
   * r &=& A / s \\
   * A &=& 1/4 * \sqrt{(a + b + c) * (a - b + c) * (a + b - c) (-a + b + c)} \\
   * s &=& \frac{a + b + c}{2}
   * @f}
   */

  Real a, b, c;
  a = distance_2d(coord1, coord2);
  b = distance_2d(coord2, coord3);
  c = distance_2d(coord1, coord3);

  Real s;
  s = (a + b + c) * 0.5;

  return sqrt((s - a) * (s - b) * (s - c) / s);
}

/* -------------------------------------------------------------------------- */
inline Real Math::distance_3d(const Real * x, const Real * y) {
  return sqrt((y[0] - x[0]) * (y[0] - x[0]) + (y[1] - x[1]) * (y[1] - x[1]) +
              (y[2] - x[2]) * (y[2] - x[2]));
}

/* -------------------------------------------------------------------------- */
inline Real Math::tetrahedron_volume(const Real * coord1, const Real * coord2,
                                     const Real * coord3, const Real * coord4) {
  Real xx[9], vol;

  xx[0] = coord2[0];
  xx[1] = coord2[1];
  xx[2] = coord2[2];
  xx[3] = coord3[0];
  xx[4] = coord3[1];
  xx[5] = coord3[2];
  xx[6] = coord4[0];
  xx[7] = coord4[1];
  xx[8] = coord4[2];
  vol = det3(xx);

  xx[0] = coord1[0];
  xx[1] = coord1[1];
  xx[2] = coord1[2];
  xx[3] = coord3[0];
  xx[4] = coord3[1];
  xx[5] = coord3[2];
  xx[6] = coord4[0];
  xx[7] = coord4[1];
  xx[8] = coord4[2];
  vol -= det3(xx);

  xx[0] = coord1[0];
  xx[1] = coord1[1];
  xx[2] = coord1[2];
  xx[3] = coord2[0];
  xx[4] = coord2[1];
  xx[5] = coord2[2];
  xx[6] = coord4[0];
  xx[7] = coord4[1];
  xx[8] = coord4[2];
  vol += det3(xx);

  xx[0] = coord1[0];
  xx[1] = coord1[1];
  xx[2] = coord1[2];
  xx[3] = coord2[0];
  xx[4] = coord2[1];
  xx[5] = coord2[2];
  xx[6] = coord3[0];
  xx[7] = coord3[1];
  xx[8] = coord3[2];
  vol -= det3(xx);

  vol /= 6;

  return vol;
}

/* -------------------------------------------------------------------------- */
inline Real Math::tetrahedron_inradius(const Real * coord1, const Real * coord2,
                                       const Real * coord3,
                                       const Real * coord4) {

  Real l12, l13, l14, l23, l24, l34;
  l12 = distance_3d(coord1, coord2);
  l13 = distance_3d(coord1, coord3);
  l14 = distance_3d(coord1, coord4);
  l23 = distance_3d(coord2, coord3);
  l24 = distance_3d(coord2, coord4);
  l34 = distance_3d(coord3, coord4);

  Real s1, s2, s3, s4;

  s1 = (l12 + l23 + l13) * 0.5;
  s1 = sqrt(s1 * (s1 - l12) * (s1 - l23) * (s1 - l13));

  s2 = (l12 + l24 + l14) * 0.5;
  s2 = sqrt(s2 * (s2 - l12) * (s2 - l24) * (s2 - l14));

  s3 = (l23 + l34 + l24) * 0.5;
  s3 = sqrt(s3 * (s3 - l23) * (s3 - l34) * (s3 - l24));

  s4 = (l13 + l34 + l14) * 0.5;
  s4 = sqrt(s4 * (s4 - l13) * (s4 - l34) * (s4 - l14));

  Real volume = Math::tetrahedron_volume(coord1, coord2, coord3, coord4);

  return 3 * volume / (s1 + s2 + s3 + s4);
}

/* -------------------------------------------------------------------------- */
inline void Math::barycenter(const Real * coord, UInt nb_points,
                             UInt spatial_dimension, Real * barycenter) {
  memset(barycenter, 0, spatial_dimension * sizeof(Real));
  for (UInt n = 0; n < nb_points; ++n) {
    UInt offset = n * spatial_dimension;
    for (UInt i = 0; i < spatial_dimension; ++i) {
      barycenter[i] += coord[offset + i] / (Real)nb_points;
    }
  }
}

/* -------------------------------------------------------------------------- */
inline void Math::vector_2d(const Real * x, const Real * y, Real * res) {
  res[0] = y[0] - x[0];
  res[1] = y[1] - x[1];
}

/* -------------------------------------------------------------------------- */
inline void Math::vector_3d(const Real * x, const Real * y, Real * res) {
  res[0] = y[0] - x[0];
  res[1] = y[1] - x[1];
  res[2] = y[2] - x[2];
}

/* -------------------------------------------------------------------------- */
/// Combined absolute and relative tolerance test proposed in
/// Real-time collision detection by C. Ericson (2004)
inline bool Math::are_float_equal(const Real x, const Real y) {
  Real abs_max = std::max(std::abs(x), std::abs(y));
  abs_max = std::max(abs_max, Real(1.));
  return std::abs(x - y) <= (tolerance * abs_max);
}

/* -------------------------------------------------------------------------- */
inline bool Math::isnan(Real x) {
#if defined(__INTEL_COMPILER)
#pragma warning(push)
#pragma warning(disable : 1572)
#endif // defined(__INTEL_COMPILER)

  // x = x return false means x = quiet_NaN
  return !(x == x);

#if defined(__INTEL_COMPILER)
#pragma warning(pop)
#endif // defined(__INTEL_COMPILER)
}

/* -------------------------------------------------------------------------- */
inline bool Math::are_vector_equal(UInt n, Real * x, Real * y) {
  bool test = true;
  for (UInt i = 0; i < n; ++i) {
    test &= are_float_equal(x[i], y[i]);
  }

  return test;
}

/* -------------------------------------------------------------------------- */
inline bool Math::intersects(Real x_min, Real x_max, Real y_min, Real y_max) {
  return not((x_max < y_min) or (x_min > y_max));
}

/* -------------------------------------------------------------------------- */
inline bool Math::is_in_range(Real a, Real x_min, Real x_max) {
  return ((a >= x_min) and (a <= x_max));
}

/* -------------------------------------------------------------------------- */
template <UInt p, typename T> inline T Math::pow(T x) {
  return (pow<p - 1, T>(x) * x);
}
template <> inline UInt Math::pow<0, UInt>(__attribute__((unused)) UInt x) {
  return (1);
}
template <> inline Real Math::pow<0, Real>(__attribute__((unused)) Real x) {
  return (1.);
}

/* -------------------------------------------------------------------------- */

template <class Functor>
Real Math::NewtonRaphson::solve(const Functor & funct, Real x_0) {
  Real x = x_0;
  Real f_x = funct.f(x);
  UInt iter = 0;
  while (std::abs(f_x) > this->tolerance && iter < this->max_iteration) {
    x -= f_x / funct.f_prime(x);
    f_x = funct.f(x);
    iter++;
  }

  AKANTU_DEBUG_ASSERT(iter < this->max_iteration,
                      "Newton Raphson ("
                          << funct.name << ") solve did not converge in "
                          << this->max_iteration << " iterations (tolerance: "
                          << this->tolerance << ")");

  return x;
}
} // namespace akantu
