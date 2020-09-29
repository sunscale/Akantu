/**
 * @file   aka_math.hh
 *
 * @author Ramin Aghababaei <ramin.aghababaei@epfl.ch>
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
 * @author Peter Spijker <peter.spijker@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Aug 04 2010
 * @date last modification: Mon Sep 11 2017
 *
 * @brief  mathematical operations
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
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */
#include <utility>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_AKA_MATH_H_
#define AKANTU_AKA_MATH_H_

namespace akantu {
/* -------------------------------------------------------------------------- */

namespace Math {
  /// tolerance for functions that need one
  extern Real tolerance; // NOLINT

  /* ------------------------------------------------------------------------ */
  /* Matrix algebra                                                           */
  /* ------------------------------------------------------------------------ */
  /// @f$ y = A*x @f$
  void matrix_vector(UInt m, UInt n, const Array<Real> & A,
                     const Array<Real> & x, Array<Real> & y, Real alpha = 1.);

  /// @f$ y = A*x @f$
  inline void matrix_vector(UInt m, UInt n, Real * A, Real * x, Real * y,
                            Real alpha = 1.);

  /// @f$ y = A^t*x @f$
  inline void matrixt_vector(UInt m, UInt n, Real * A, Real * x, Real * y,
                             Real alpha = 1.);

  /// @f$ C = A*B @f$
  void matrix_matrix(UInt m, UInt n, UInt k, const Array<Real> & A,
                     const Array<Real> & B, Array<Real> & C, Real alpha = 1.);

  /// @f$ C = A*B^t @f$
  void matrix_matrixt(UInt m, UInt n, UInt k, const Array<Real> & A,
                      const Array<Real> & B, Array<Real> & C, Real alpha = 1.);

  /// @f$ C = A*B @f$
  inline void matrix_matrix(UInt m, UInt n, UInt k, Real * A, Real * B,
                            Real * C, Real alpha = 1.);

  /// @f$ C = A^t*B @f$
  inline void matrixt_matrix(UInt m, UInt n, UInt k, Real * A, Real * B,
                             Real * C, Real alpha = 1.);

  /// @f$ C = A*B^t @f$
  inline void matrix_matrixt(UInt m, UInt n, UInt k, Real * A, Real * B,
                             Real * C, Real alpha = 1.);

  /// @f$ C = A^t*B^t @f$
  inline void matrixt_matrixt(UInt m, UInt n, UInt k, Real * A, Real * B,
                              Real * C, Real alpha = 1.);

  template <bool tr_A, bool tr_B>
  inline void matMul(UInt m, UInt n, UInt k, Real alpha, Real * A, Real * B,
                     Real beta, Real * C);

  template <bool tr_A>
  inline void matVectMul(UInt m, UInt n, Real alpha, Real * A, Real * x,
                         Real beta, Real * y);

  inline void aXplusY(UInt n, Real alpha, Real * x, Real * y);

  inline void matrix33_eigenvalues(Real * A, Real * Adiag);

  inline void matrix22_eigenvalues(Real * A, Real * Adiag);
  template <UInt dim> inline void eigenvalues(Real * A, Real * d);

  /// solve @f$ A x = \Lambda x @f$ and return d and V such as @f$ A V[i:] =
  /// d[i] V[i:]@f$
  template <typename T> void matrixEig(UInt n, T * A, T * d, T * V = nullptr);

  /// determinent of a 2x2 matrix
  Real det2(const Real * mat);
  /// determinent of a 3x3 matrix
  Real det3(const Real * mat);
  /// determinent of a nxn matrix
  template <UInt n> Real det(const Real * mat);
  /// determinent of a nxn matrix
  template <typename T> T det(UInt n, const T * A);

  /// inverse a nxn matrix
  template <UInt n> inline void inv(const Real * A, Real * inv);
  /// inverse a nxn matrix
  template <typename T> inline void inv(UInt n, const T * A, T * inv);
  /// inverse a 3x3 matrix
  inline void inv3(const Real * mat, Real * inv);
  /// inverse a 2x2 matrix
  inline void inv2(const Real * mat, Real * inv);

  /// solve A x = b using a LU factorization
  template <typename T>
  inline void solve(UInt n, const T * A, T * x, const T * b);

  /// return the double dot product between 2 tensors in 2d
  inline Real matrixDoubleDot22(Real * A, Real * B);

  /// return the double dot product between 2 tensors in 3d
  inline Real matrixDoubleDot33(Real * A, Real * B);

  /// extension of the double dot product to two 2nd order tensor in dimension n
  inline Real matrixDoubleDot(UInt n, Real * A, Real * B);

  /* ------------------------------------------------------------------------ */
  /* Array algebra                                                            */
  /* ------------------------------------------------------------------------ */
  /// vector cross product
  inline void vectorProduct3(const Real * v1, const Real * v2, Real * res);

  /// normalize a vector
  inline void normalize2(Real * v);

  /// normalize a vector
  inline void normalize3(Real * v);

  /// return norm of a 2-vector
  inline Real norm2(const Real * v);

  /// return norm of a 3-vector
  inline Real norm3(const Real * v);

  /// return norm of a vector
  inline Real norm(UInt n, const Real * v);

  /// return the dot product between 2 vectors in 2d
  inline Real vectorDot2(const Real * v1, const Real * v2);

  /// return the dot product between 2 vectors in 3d
  inline Real vectorDot3(const Real * v1, const Real * v2);

  /// return the dot product between 2 vectors
  inline Real vectorDot(Real * v1, Real * v2, UInt n);

  /* ------------------------------------------------------------------------ */
  /* Geometry                                                                 */
  /* ------------------------------------------------------------------------ */
  /// compute normal a normal to a vector
  inline void normal2(const Real * vec, Real * normal);

  /// compute normal a normal to a vector
  inline void normal3(const Real * vec1, const Real * vec2, Real * normal);

  /// compute the tangents to an array of normal vectors
  void compute_tangents(const Array<Real> & normals, Array<Real> & tangents);

  /// distance in 2D between x and y
  inline Real distance_2d(const Real * x, const Real * y);

  /// distance in 3D between x and y
  inline Real distance_3d(const Real * x, const Real * y);

  /// radius of the in-circle of a triangle
  inline Real triangle_inradius(const Real * coord1, const Real * coord2,
                                const Real * coord3);

  /// radius of the in-circle of a tetrahedron
  inline Real tetrahedron_inradius(const Real * coord1, const Real * coord2,
                                   const Real * coord3, const Real * coord4);

  /// volume of a tetrahedron
  inline Real tetrahedron_volume(const Real * coord1, const Real * coord2,
                                 const Real * coord3, const Real * coord4);

  /// compute the barycenter of n points
  inline void barycenter(const Real * coord, UInt nb_points,
                         UInt spatial_dimension, Real * barycenter);

  /// vector between x and y
  inline void vector_2d(const Real * x, const Real * y, Real * res);

  /// vector pointing from x to y in 3 spatial dimension
  inline void vector_3d(const Real * x, const Real * y, Real * res);

  /// test if two scalar are equal within a given tolerance
  inline bool are_float_equal(Real x, Real y);

  /// test if two vectors are equal within a given tolerance
  inline bool are_vector_equal(UInt n, Real * x, Real * y);

#ifdef isnan
#error                                                                         \
    "You probably  included <math.h> which  is incompatible with aka_math  please use\
<cmath> or add a \"#undef isnan\" before akantu includes"
#endif
  /// test if a real is a NaN
  inline bool isnan(Real x);

  /// test if the line x and y intersects each other
  inline bool intersects(Real x_min, Real x_max, Real y_min, Real y_max);

  /// test if a is in the range [x_min, x_max]
  inline bool is_in_range(Real a, Real x_min, Real x_max);

  inline Real getTolerance() { return Math::tolerance; };
  inline void setTolerance(Real tol) { Math::tolerance = tol; };

  template <UInt p, typename T> inline T pow(T x);

  template <class T1, class T2,
            std::enable_if_t<std::is_integral<T1>::value and
                             std::is_integral<T2>::value> * = nullptr>
  inline Real kronecker(T1 i, T2 j) {
    return static_cast<Real>(i == j);
  }

  /// reduce all the values of an array, the summation is done in place and the
  /// array is modified
  Real reduce(Array<Real> & array);

  class NewtonRaphson {
  public:
    NewtonRaphson(Real tolerance, Real max_iteration)
        : tolerance(tolerance), max_iteration(max_iteration) {}

    template <class Functor> Real solve(const Functor & funct, Real x_0);

  private:
    Real tolerance;
    Real max_iteration;
  };

  struct NewtonRaphsonFunctor {
    explicit NewtonRaphsonFunctor(std::string name) : name(std::move(name)) {}

    virtual ~NewtonRaphsonFunctor() = default;

    NewtonRaphsonFunctor(const NewtonRaphsonFunctor & other) = default;
    NewtonRaphsonFunctor(NewtonRaphsonFunctor && other) noexcept = default;
    NewtonRaphsonFunctor &
    operator=(const NewtonRaphsonFunctor & other) = default;
    NewtonRaphsonFunctor &
    operator=(NewtonRaphsonFunctor && other) noexcept = default;

    virtual Real f(Real x) const = 0;
    virtual Real f_prime(Real x) const = 0;
    std::string name;
  };
} // namespace Math
} // namespace akantu
/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "aka_math_tmpl.hh"

#endif /* AKANTU_AKA_MATH_H_ */
