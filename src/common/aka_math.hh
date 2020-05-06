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

#ifndef __AKANTU_AKA_MATH_H__
#define __AKANTU_AKA_MATH_H__

namespace akantu {
/* -------------------------------------------------------------------------- */


class Math {
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Matrix algebra                                                           */
  /* ------------------------------------------------------------------------ */
  /// @f$ y = A*x @f$
  static void matrix_vector(UInt m, UInt n, const Array<Real> & A,
                            const Array<Real> & x, Array<Real> & y,
                            Real alpha = 1.);

  /// @f$ y = A*x @f$
  static inline void matrix_vector(UInt m, UInt n, Real * A, Real * x, Real * y,
                                   Real alpha = 1.);

  /// @f$ y = A^t*x @f$
  static inline void matrixt_vector(UInt m, UInt n, Real * A, Real * x,
                                    Real * y, Real alpha = 1.);

  /// @f$ C = A*B @f$
  static void matrix_matrix(UInt m, UInt n, UInt k, const Array<Real> & A,
                            const Array<Real> & B, Array<Real> & C,
                            Real alpha = 1.);

  /// @f$ C = A*B^t @f$
  static void matrix_matrixt(UInt m, UInt n, UInt k,
                             const Array<Real> & A,
                             const Array<Real> & B, Array<Real> & C,
                             Real alpha = 1.);

  /// @f$ C = A*B @f$
  static inline void matrix_matrix(UInt m, UInt n, UInt k, Real * A, Real * B,
                                   Real * C, Real alpha = 1.);

  /// @f$ C = A^t*B @f$
  static inline void matrixt_matrix(UInt m, UInt n, UInt k, Real * A, Real * B,
                                    Real * C, Real alpha = 1.);

  /// @f$ C = A*B^t @f$
  static inline void matrix_matrixt(UInt m, UInt n, UInt k, Real * A, Real * B,
                                    Real * C, Real alpha = 1.);

  /// @f$ C = A^t*B^t @f$
  static inline void matrixt_matrixt(UInt m, UInt n, UInt k, Real * A, Real * B,
                                     Real * C, Real alpha = 1.);

  template <bool tr_A, bool tr_B>
  static inline void matMul(UInt m, UInt n, UInt k, Real alpha, Real * A,
                            Real * B, Real beta, Real * C);

  template <bool tr_A>
  static inline void matVectMul(UInt m, UInt n, Real alpha, Real * A, Real * x,
                                Real beta, Real * y);

  static inline void aXplusY(UInt n, Real alpha, Real * x, Real * y);

  static inline void matrix33_eigenvalues(Real * A, Real * Adiag);

  static inline void matrix22_eigenvalues(Real * A, Real * Adiag);
  template <UInt dim> static inline void eigenvalues(Real * A, Real * d);

  /// solve @f$ A x = \Lambda x @f$ and return d and V such as @f$ A V[i:] =
  /// d[i] V[i:]@f$
  template <typename T>
  static void matrixEig(UInt n, T * A, T * d, T * V = nullptr);

  /// determinent of a 2x2 matrix
  static inline Real det2(const Real * mat);
  /// determinent of a 3x3 matrix
  static inline Real det3(const Real * mat);
  /// determinent of a nxn matrix
  template <UInt n> static inline Real det(const Real * mat);
  /// determinent of a nxn matrix
  template <typename T> static inline T det(UInt n, const T * mat);

  /// inverse a nxn matrix
  template <UInt n> static inline void inv(const Real * mat, Real * inv);
  /// inverse a nxn matrix
  template <typename T> static inline void inv(UInt n, const T * mat, T * inv);
  /// inverse a 3x3 matrix
  static inline void inv3(const Real * mat, Real * inv);
  /// inverse a 2x2 matrix
  static inline void inv2(const Real * mat, Real * inv);

  /// solve A x = b using a LU factorization
  template <typename T>
  static inline void solve(UInt n, const T * A, T * x, const T * b);

  /// return the double dot product between 2 tensors in 2d
  static inline Real matrixDoubleDot22(Real * A, Real * B);

  /// return the double dot product between 2 tensors in 3d
  static inline Real matrixDoubleDot33(Real * A, Real * B);

  /// extension of the double dot product to two 2nd order tensor in dimension n
  static inline Real matrixDoubleDot(UInt n, Real * A, Real * B);

  /* ------------------------------------------------------------------------ */
  /* Array algebra                                                            */
  /* ------------------------------------------------------------------------ */
  /// vector cross product
  static inline void vectorProduct3(const Real * v1, const Real * v2,
                                    Real * res);

  /// normalize a vector
  static inline void normalize2(Real * v);

  /// normalize a vector
  static inline void normalize3(Real * v);

  /// return norm of a 2-vector
  static inline Real norm2(const Real * v);

  /// return norm of a 3-vector
  static inline Real norm3(const Real * v);

  /// return norm of a vector
  static inline Real norm(UInt n, const Real * v);

  /// return the dot product between 2 vectors in 2d
  static inline Real vectorDot2(const Real * v1, const Real * v2);

  /// return the dot product between 2 vectors in 3d
  static inline Real vectorDot3(const Real * v1, const Real * v2);

  /// return the dot product between 2 vectors
  static inline Real vectorDot(Real * v1, Real * v2, UInt n);

  /* ------------------------------------------------------------------------ */
  /* Geometry                                                                 */
  /* ------------------------------------------------------------------------ */
  /// compute normal a normal to a vector
  static inline void normal2(const Real * v1, Real * res);

  /// compute normal a normal to a vector
  static inline void normal3(const Real * v1, const Real * v2, Real * res);

  /// compute the tangents to an array of normal vectors
  static void compute_tangents(const Array<Real> & normals,
                               Array<Real> & tangents);

  /// distance in 2D between x and y
  static inline Real distance_2d(const Real * x, const Real * y);

  /// distance in 3D between x and y
  static inline Real distance_3d(const Real * x, const Real * y);

  /// radius of the in-circle of a triangle
  static inline Real triangle_inradius(const Real * coord1, const Real * coord2,
                                       const Real * coord3);

  /// radius of the in-circle of a tetrahedron
  static inline Real tetrahedron_inradius(const Real * coord1,
                                          const Real * coord2,
                                          const Real * coord3,
                                          const Real * coord4);

  /// volume of a tetrahedron
  static inline Real tetrahedron_volume(const Real * coord1,
                                        const Real * coord2,
                                        const Real * coord3,
                                        const Real * coord4);

  /// compute the barycenter of n points
  static inline void barycenter(const Real * coord, UInt nb_points,
                                UInt spatial_dimension, Real * barycenter);

  /// vector between x and y
  static inline void vector_2d(const Real * x, const Real * y, Real * vec);

  /// vector pointing from x to y in 3 spatial dimension
  static inline void vector_3d(const Real * x, const Real * y, Real * vec);

  /// test if two scalar are equal within a given tolerance
  static inline bool are_float_equal(Real x, Real y);

  /// test if two vectors are equal within a given tolerance
  static inline bool are_vector_equal(UInt n, Real * x, Real * y);

#ifdef isnan
#error                                                                         \
    "You probably  included <math.h> which  is incompatible with aka_math  please use\
<cmath> or add a \"#undef isnan\" before akantu includes"
#endif
  /// test if a real is a NaN
  static inline bool isnan(Real x);

  /// test if the line x and y intersects each other
  static inline bool intersects(Real x_min, Real x_max, Real y_min, Real y_max);

  /// test if a is in the range [x_min, x_max]
  static inline bool is_in_range(Real a, Real x_min, Real x_max);

  static inline Real getTolerance() { return tolerance; };
  static inline void setTolerance(Real tol) { tolerance = tol; };

  template <UInt p, typename T> static inline T pow(T x);

  /// reduce all the values of an array, the summation is done in place and the
  /// array is modified
  static Real reduce(Array<Real> & array);

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
    virtual Real f(Real x) const = 0;
    virtual Real f_prime(Real x) const = 0;
    std::string name;
  };

private:
  /// tolerance for functions that need one
  static Real tolerance;
};

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "aka_math_tmpl.hh"

#endif /* __AKANTU_AKA_MATH_H__ */
