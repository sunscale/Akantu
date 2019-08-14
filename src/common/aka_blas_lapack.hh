/**
 * @file   aka_blas_lapack.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Mar 06 2013
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Interface of the Fortran BLAS/LAPACK libraries
 *
 * @section LICENSE
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

#ifndef __AKANTU_AKA_BLAS_LAPACK_HH__
#define __AKANTU_AKA_BLAS_LAPACK_HH__

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_BLAS
#include "aka_fortran_mangling.hh"
extern "C" {

/* ------------------------------------------------------------------------ */
/* Double precision                                                         */
/* ------------------------------------------------------------------------ */
// LEVEL 1
double AKA_FC_GLOBAL(ddot, DDOT)(int *, double *, int *, double *, int *);
void AKA_FC_GLOBAL(daxpy, DAXPY)(int *, double *, double *, int *, double *,
                                 int *);

// LEVEL 2
void AKA_FC_GLOBAL(dgemv, DGEMV)(char *, int *, int *, double *, double *,
                                 int *, double *, int *, double *, double *,
                                 int *);
// LEVEL 3
void AKA_FC_GLOBAL(dgemm, DGEMM)(char *, char *, int *, int *, int *, double *,
                                 double *, int *, double *, int *, double *,
                                 double *, int *);
/* ------------------------------------------------------------------------ */
/* Simple precision                                                         */
/* ------------------------------------------------------------------------ */
// LEVEL 1
float AKA_FC_GLOBAL(sdot, SDOT)(int *, float *, int *, float *, int *);
void AKA_FC_GLOBAL(saxpy, SAXPY)(int *, float *, float *, int *, float *,
                                 int *);

// LEVEL 2
void AKA_FC_GLOBAL(sgemv, SGEMV)(char *, int *, int *, float *, float *, int *,
                                 float *, int *, float *, float *, int *);
// LEVEL 3
void AKA_FC_GLOBAL(sgemm, SGEMM)(char *, char *, int *, int *, int *, float *,
                                 float *, int *, float *, int *, float *,
                                 float *, int *);
}
#endif

namespace akantu {

#define AKANTU_WARNING_IGNORE_UNUSED_PARAMETER
#include "aka_warning.hh"

/// Wrapper around the S/DDOT BLAS function that returns the dot product of two
/// vectors
template <typename T>
inline T aka_dot(int * n, T * x, int * incx, T * y, int * incy) {
  AKANTU_ERROR(debug::demangle(typeid(T).name())
               << "is not a type recognized, or you didn't activated "
                  "BLAS in the compilation options!");
}

/// Wrapper around the S/DAXPY BLAS function that computes \f$y := \alpha x +
/// y\f$
template <typename T>
inline void aka_axpy(int * n, T * alpha, T * x, int * incx, T * y, int * incy) {
  AKANTU_ERROR(debug::demangle(typeid(T).name())
               << "is not a type recognized, or you didn't activated "
                  "BLAS in the compilation options!");
}

/// Wrapper around the S/DGEMV BLAS function that computes matrix-vector product
/// \f$y := \alpha A^{(T)}x + \beta y \f$
template <typename T>
inline void aka_gemv(char * trans, int * m, int * n, T * alpha, T * a,
                     int * lda, T * x, int * incx, T * beta, T * y,
                     int * incy) {
  AKANTU_ERROR(debug::demangle(typeid(T).name())
               << "is not a type recognized, or you didn't activated "
                  "BLAS in the compilation options!");
}

/// Wrapper around the S/DGEMM BLAS function that computes the product of two
/// matrices \f$C := \alpha A^{(T)} B^{(T)} + \beta C \f$
template <typename T>
inline void aka_gemm(char * transa, char * transb, int * m, int * n, int * k,
                     T * alpha, T * a, int * lda, T * b, int * ldb, T * beta,
                     T * c, int * ldc) {
  AKANTU_ERROR(debug::demangle(typeid(T).name())
               << "is not a type recognized, or you didn't activated "
                  "BLAS in the compilation options!");
}

#if defined(AKANTU_USE_BLAS)
template <>
inline double aka_dot<double>(int * n, double * x, int * incx, double * y,
                              int * incy) {
  return AKA_FC_GLOBAL(ddot, DDOT)(n, x, incx, y, incy);
}

template <>
inline void aka_axpy(int * n, double * alpha, double * x, int * incx,
                     double * y, int * incy) {
  return AKA_FC_GLOBAL(daxpy, DAXPY)(n, alpha, x, incx, y, incy);
}

template <>
inline void aka_gemv<double>(char * trans, int * m, int * n, double * alpha,
                             double * a, int * lda, double * x, int * incx,
                             double * beta, double * y, int * incy) {
  return AKA_FC_GLOBAL(dgemv, DGEMV)(trans, m, n, alpha, a, lda, x, incx, beta,
                                     y, incy);
}

template <>
inline void aka_gemm<double>(char * transa, char * transb, int * m, int * n,
                             int * k, double * alpha, double * a, int * lda,
                             double * b, int * ldb, double * beta, double * c,
                             int * ldc) {
  AKA_FC_GLOBAL(dgemm, DGEMM)
  (transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

template <>
inline float aka_dot<float>(int * n, float * x, int * incx, float * y,
                            int * incy) {
  return AKA_FC_GLOBAL(sdot, SDOT)(n, x, incx, y, incy);
}

template <>
inline void aka_axpy(int * n, float * alpha, float * x, int * incx, float * y,
                     int * incy) {
  return AKA_FC_GLOBAL(daxpy, DAXPY)(n, alpha, x, incx, y, incy);
}

template <>
inline void aka_gemv<float>(char * trans, int * m, int * n, float * alpha,
                            float * a, int * lda, float * x, int * incx,
                            float * beta, float * y, int * incy) {
  AKA_FC_GLOBAL(sgemv, SGEMV)
  (trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
}

template <>
inline void aka_gemm<float>(char * transa, char * transb, int * m, int * n,
                            int * k, float * alpha, float * a, int * lda,
                            float * b, int * ldb, float * beta, float * c,
                            int * ldc) {
  AKA_FC_GLOBAL(sgemm, SGEMM)
  (transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

#endif

} // namespace akantu

#ifdef AKANTU_USE_LAPACK
#include "aka_fortran_mangling.hh"
extern "C" {
/* ------------------------------------------------------------------------ */
/* Double general matrix                                                    */
/* ------------------------------------------------------------------------ */
/// compute the eigenvalues/vectors
void AKA_FC_GLOBAL(dgeev, DGEEV)(char * jobvl, char * jobvr, int * n,
                                 double * a, int * lda, double * wr,
                                 double * wi, double * vl, int * ldvl,
                                 double * vr, int * ldvr, double * work,
                                 int * lwork, int * info);

/// LU decomposition of a general matrix
void AKA_FC_GLOBAL(dgetrf, DGETRF)(int * m, int * n, double * a, int * lda,
                                   int * ipiv, int * info);

/// generate inverse of a matrix given its LU decomposition
void AKA_FC_GLOBAL(dgetri, DGETRI)(int * n, double * a, int * lda, int * ipiv,
                                   double * work, int * lwork, int * info);

/// solving A x = b using a LU factorization
void AKA_FC_GLOBAL(dgetrs, DGETRS)(char * trans, int * n, int * nrhs,
                                   double * A, int * lda, int * ipiv,
                                   double * b, int * ldb, int * info);

/* ------------------------------------------------------------------------ */
/* Simple general matrix                                                    */
/* ------------------------------------------------------------------------ */
/// compute the eigenvalues/vectors
void AKA_FC_GLOBAL(sgeev, SGEEV)(char * jobvl, char * jobvr, int * n, float * a,
                                 int * lda, float * wr, float * wi, float * vl,
                                 int * ldvl, float * vr, int * ldvr,
                                 float * work, int * lwork, int * info);

/// LU decomposition of a general matrix
void AKA_FC_GLOBAL(sgetrf, SGETRF)(int * m, int * n, float * a, int * lda,
                                   int * ipiv, int * info);

/// generate inverse of a matrix given its LU decomposition
void AKA_FC_GLOBAL(sgetri, SGETRI)(int * n, float * a, int * lda, int * ipiv,
                                   float * work, int * lwork, int * info);

/// solving A x = b using a LU factorization
void AKA_FC_GLOBAL(sgetrs, SGETRS)(char * trans, int * n, int * nrhs, float * A,
                                   int * lda, int * ipiv, float * b, int * ldb,
                                   int * info);
}
#endif // AKANTU_USE_LAPACK

namespace akantu {

/// Wrapper around the S/DGEEV BLAS function that computes the eigenvalues and
/// eigenvectors of a matrix
template <typename T>
inline void aka_geev(char * jobvl, char * jobvr, int * n, T * a, int * lda,
                     T * wr, T * wi, T * vl, int * ldvl, T * vr, int * ldvr,
                     T * work, int * lwork, int * info) {
  AKANTU_ERROR(debug::demangle(typeid(T).name())
               << "is not a type recognized, or you didn't activated "
                  "LAPACK in the compilation options!");
}

/// Wrapper around the S/DGETRF BLAS function that computes the LU decomposition
/// of a matrix
template <typename T>
inline void aka_getrf(int * m, int * n, T * a, int * lda, int * ipiv,
                      int * info) {
  AKANTU_ERROR(debug::demangle(typeid(T).name())
               << "is not a type recognized, or you didn't activated "
                  "LAPACK in the compilation options!");
}

/// Wrapper around the S/DGETRI BLAS function that computes the inverse of a
/// matrix given its LU decomposition
template <typename T>
inline void aka_getri(int * n, T * a, int * lda, int * ipiv, T * work,
                      int * lwork, int * info) {
  AKANTU_ERROR(debug::demangle(typeid(T).name())
               << "is not a type recognized, or you didn't activated "
                  "LAPACK in the compilation options!");
}

/// Wrapper around the S/DGETRS BLAS function that solves \f$A^{(T)}x = b\f$
/// using LU decomposition
template <typename T>
inline void aka_getrs(char * trans, int * n, int * nrhs, T * A, int * lda,
                      int * ipiv, T * b, int * ldb, int * info) {
  AKANTU_ERROR(debug::demangle(typeid(T).name())
               << "is not a type recognized, or you didn't activated "
                  "LAPACK in the compilation options!");
}

#include "aka_warning_restore.hh"

#ifdef AKANTU_USE_LAPACK
template <>
inline void aka_geev<double>(char * jobvl, char * jobvr, int * n, double * a,
                             int * lda, double * wr, double * wi, double * vl,
                             int * ldvl, double * vr, int * ldvr, double * work,
                             int * lwork, int * info) {
  AKA_FC_GLOBAL(dgeev, DGEEV)
  (jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info);
}

template <>
inline void aka_getrf<double>(int * m, int * n, double * a, int * lda,
                              int * ipiv, int * info) {
  AKA_FC_GLOBAL(dgetrf, DGETRF)(m, n, a, lda, ipiv, info);
}

template <>
inline void aka_getri<double>(int * n, double * a, int * lda, int * ipiv,
                              double * work, int * lwork, int * info) {
  AKA_FC_GLOBAL(dgetri, DGETRI)(n, a, lda, ipiv, work, lwork, info);
}

template <>
inline void aka_getrs<double>(char * trans, int * n, int * nrhs, double * A,
                              int * lda, int * ipiv, double * b, int * ldb,
                              int * info) {
  AKA_FC_GLOBAL(dgetrs, DGETRS)(trans, n, nrhs, A, lda, ipiv, b, ldb, info);
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
template <>
inline void aka_geev<float>(char * jobvl, char * jobvr, int * n, float * a,
                            int * lda, float * wr, float * wi, float * vl,
                            int * ldvl, float * vr, int * ldvr, float * work,
                            int * lwork, int * info) {
  AKA_FC_GLOBAL(sgeev, SGEEV)
  (jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info);
}

template <>
inline void aka_getrf<float>(int * m, int * n, float * a, int * lda, int * ipiv,
                             int * info) {
  AKA_FC_GLOBAL(sgetrf, SGETRF)(m, n, a, lda, ipiv, info);
}

template <>
inline void aka_getri<float>(int * n, float * a, int * lda, int * ipiv,
                             float * work, int * lwork, int * info) {
  AKA_FC_GLOBAL(sgetri, SGETRI)(n, a, lda, ipiv, work, lwork, info);
}

template <>
inline void aka_getrs<float>(char * trans, int * n, int * nrhs, float * A,
                             int * lda, int * ipiv, float * b, int * ldb,
                             int * info) {
  AKA_FC_GLOBAL(sgetrs, SGETRS)(trans, n, nrhs, A, lda, ipiv, b, ldb, info);
}
#endif

} // namespace akantu

#endif /* __AKANTU_AKA_BLAS_LAPACK_HH__ */
