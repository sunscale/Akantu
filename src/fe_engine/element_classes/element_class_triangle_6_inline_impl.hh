/**
 * @file   element_class_triangle_6_inline_impl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jul 16 2010
 * @date last modification: Wed Oct 11 2017
 *
 * @brief  Specialization of the element_class class for the type _triangle_6
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 terms  of the  GNU Lesser  General Public  License as published by  the Free
 Software Foundation, either version 3 of the License, or (at your option) any
 later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * @verbatim
        \eta
          ^
          |
          x 2
          | `
          |   `
          |  .  `
          |   q2  `
        5 x          x 4
          |           `
          |             `
          |  .q0     q1.  `
          |                 `
          x---------x---------x-----> \xi
          0         3         1
 @endverbatim
 *
 *
 * @f[
 * \begin{array}{ll}
 *   \xi_{0}  = 0   &  \eta_{0} = 0   \\
 *   \xi_{1}  = 1   &  \eta_{1} = 0   \\
 *   \xi_{2}  = 0   &  \eta_{2} = 1   \\
 *   \xi_{3}  = 1/2 &  \eta_{3} = 0   \\
 *   \xi_{4}  = 1/2 &  \eta_{4} = 1/2 \\
 *   \xi_{5}  = 0   &  \eta_{5} = 1/2
 * \end{array}
 * @f]
 *
 * @f[
 * \begin{array}{lll}
 * N1 = -(1 - \xi - \eta) (1 - 2 (1 - \xi - \eta))
 *       & \frac{\partial N1}{\partial \xi}  = 1 - 4(1 - \xi - \eta)
 *       & \frac{\partial N1}{\partial \eta} = 1 - 4(1 - \xi - \eta) \\
 * N2 = - \xi (1 - 2 \xi)
 *       & \frac{\partial N2}{\partial \xi}  = - 1 + 4 \xi
 *       & \frac{\partial N2}{\partial \eta} = 0 \\
 * N3 = - \eta (1 - 2 \eta)
 *       & \frac{\partial N3}{\partial \xi}  = 0
 *       & \frac{\partial N3}{\partial \eta} = - 1 + 4 \eta \\
 * N4 = 4 \xi (1 - \xi - \eta)
 *       & \frac{\partial N4}{\partial \xi}  = 4 (1 - 2 \xi - \eta)
 *       & \frac{\partial N4}{\partial \eta} = - 4 \xi \\
 * N5 = 4 \xi \eta
 *       & \frac{\partial N5}{\partial \xi}  = 4 \eta
 *       & \frac{\partial N5}{\partial \eta} = 4 \xi \\
 * N6 = 4 \eta (1 - \xi - \eta)
 *       & \frac{\partial N6}{\partial \xi}  = - 4 \eta
 *       & \frac{\partial N6}{\partial \eta} = 4 (1 - \xi - 2 \eta)
 * \end{array}
 * @f]
 *
 * @f{eqnarray*}{
 * \xi_{q0}  &=& 1/6 \qquad  \eta_{q0} = 1/6 \\
 * \xi_{q1}  &=& 2/3 \qquad  \eta_{q1} = 1/6 \\
 * \xi_{q2}  &=& 1/6 \qquad  \eta_{q2} = 2/3
 * @f}
 */
/* -------------------------------------------------------------------------- */
#include "element_class.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
/* -------------------------------------------------------------------------- */
AKANTU_DEFINE_ELEMENT_CLASS_PROPERTY(_triangle_6, _gt_triangle_6,
                                     _itp_lagrange_triangle_6, _ek_regular, 2,
                                     _git_triangle, 2);

/* -------------------------------------------------------------------------- */

template <>
template <class vector_type>
inline void InterpolationElement<_itp_lagrange_triangle_6>::computeShapes(
    const vector_type & natural_coords, vector_type & N) {
  /// Natural coordinates
  Real c0 =
      1 - natural_coords(0) - natural_coords(1); /// @f$ c0 = 1 - \xi - \eta @f$
  Real c1 = natural_coords(0);                   /// @f$ c1 = \xi @f$
  Real c2 = natural_coords(1);                   /// @f$ c2 = \eta @f$

  N(0) = c0 * (2 * c0 - 1.);
  N(1) = c1 * (2 * c1 - 1.);
  N(2) = c2 * (2 * c2 - 1.);
  N(3) = 4 * c0 * c1;
  N(4) = 4 * c1 * c2;
  N(5) = 4 * c2 * c0;
}
/* -------------------------------------------------------------------------- */
template <>
template <class vector_type, class matrix_type>
inline void InterpolationElement<_itp_lagrange_triangle_6>::computeDNDS(
    const vector_type & natural_coords, matrix_type & dnds) {

  /**
   * @f[
   * dnds =  \left(
   *   \begin{array}{cccccc}
   *       \frac{\partial N1}{\partial \xi}
   *     & \frac{\partial N2}{\partial \xi}
   *     & \frac{\partial N3}{\partial \xi}
   *     & \frac{\partial N4}{\partial \xi}
   *     & \frac{\partial N5}{\partial \xi}
   *     & \frac{\partial N6}{\partial \xi} \\
   *
   *       \frac{\partial N1}{\partial \eta}
   *     & \frac{\partial N2}{\partial \eta}
   *     & \frac{\partial N3}{\partial \eta}
   *     & \frac{\partial N4}{\partial \eta}
   *     & \frac{\partial N5}{\partial \eta}
   *     & \frac{\partial N6}{\partial \eta}
   *   \end{array}
   * \right)
   * @f]
   */

  /// Natural coordinates
  Real c0 =
      1 - natural_coords(0) - natural_coords(1); /// @f$ c0 = 1 - \xi - \eta @f$
  Real c1 = natural_coords(0);                   /// @f$ c1 = \xi @f$
  Real c2 = natural_coords(1);                   /// @f$ c2 = \eta @f$

  dnds(0, 0) = 1 - 4 * c0;
  dnds(0, 1) = 4 * c1 - 1.;
  dnds(0, 2) = 0.;
  dnds(0, 3) = 4 * (c0 - c1);
  dnds(0, 4) = 4 * c2;
  dnds(0, 5) = -4 * c2;

  dnds(1, 0) = 1 - 4 * c0;
  dnds(1, 1) = 0.;
  dnds(1, 2) = 4 * c2 - 1.;
  dnds(1, 3) = -4 * c1;
  dnds(1, 4) = 4 * c1;
  dnds(1, 5) = 4 * (c0 - c2);
}

/* -------------------------------------------------------------------------- */
template <>
inline void
InterpolationElement<_itp_lagrange_triangle_6>::computeSpecialJacobian(
    const Matrix<Real> & J, Real & jac) {
  Vector<Real> vprod(J.cols());
  Matrix<Real> Jt(J.transpose(), true);
  vprod.crossProduct(Jt(0), Jt(1));
  jac = vprod.norm();
}

/* -------------------------------------------------------------------------- */
template <>
inline Real
GeometricalElement<_gt_triangle_6>::getInradius(const Matrix<Real> & coord) {
  UInt triangles[4][3] = {{0, 3, 5}, {3, 1, 4}, {3, 4, 5}, {5, 4, 2}};

  Real inradius = std::numeric_limits<Real>::max();
  for (UInt t = 0; t < 4; t++) {
    Real ir = Math::triangle_inradius(coord(triangles[t][0]).storage(),
                                      coord(triangles[t][1]).storage(),
                                      coord(triangles[t][2]).storage());
    inradius = std::min(ir, inradius);
  }

  return 2. * inradius;
}

/* -------------------------------------------------------------------------- */
// template<> inline bool ElementClass<_triangle_6>::contains(const Vector<Real>
// & natural_coords) {
//   return ElementClass<_triangle_3>::contains(natural_coords);
// }
} // namespace akantu
