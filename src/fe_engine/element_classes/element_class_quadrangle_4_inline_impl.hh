/**
 * @file   element_class_quadrangle_4_inline_impl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Dec 13 2010
 * @date last modification: Wed Oct 11 2017
 *
 * @brief  Specialization of the element_class class for the type _quadrangle_4
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
 (-1,1)   |   (1,1)
     x---------x
     |    |    |
     |    |    |
   --|---------|----->  \xi
     |    |    |
     |    |    |
     x---------x
 (-1,-1)  |   (1,-1)
 @endverbatim
 *
 * @f[
 * \begin{array}{lll}
 * N1 = (1 - \xi) (1 - \eta) / 4
 *       & \frac{\partial N1}{\partial \xi}  = - (1 - \eta) / 4
 *       & \frac{\partial N1}{\partial \eta} = - (1 - \xi) / 4 \\
 * N2 = (1 + \xi) (1 - \eta) / 4 \\
 *       & \frac{\partial N2}{\partial \xi}  = (1 - \eta) / 4
 *       & \frac{\partial N2}{\partial \eta} = - (1 + \xi) / 4 \\
 * N3 = (1 + \xi) (1 + \eta) / 4 \\
 *       & \frac{\partial N3}{\partial \xi}  = (1 + \eta) / 4
 *       & \frac{\partial N3}{\partial \eta} = (1 + \xi) / 4 \\
 * N4 = (1 - \xi) (1 + \eta) / 4
 *       & \frac{\partial N4}{\partial \xi}  = - (1 + \eta) / 4
 *       & \frac{\partial N4}{\partial \eta} = (1 - \xi) / 4 \\
 * \end{array}
 * @f]
 *
 * @f{eqnarray*}{
 * \xi_{q0}  &=& 0 \qquad  \eta_{q0} = 0
 * @f}
 */
/* -------------------------------------------------------------------------- */
#include "element_class.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
/* -------------------------------------------------------------------------- */
AKANTU_DEFINE_ELEMENT_CLASS_PROPERTY(_quadrangle_4, _gt_quadrangle_4,
                                     _itp_lagrange_quadrangle_4, _ek_regular, 2,
                                     _git_segment, 2);

/* -------------------------------------------------------------------------- */
template <>
template <class vector_type>
inline void InterpolationElement<_itp_lagrange_quadrangle_4>::computeShapes(
    const vector_type & c, vector_type & N) {
  N(0) = 1. / 4. * (1. - c(0)) * (1. - c(1)); /// N1(q_0)
  N(1) = 1. / 4. * (1. + c(0)) * (1. - c(1)); /// N2(q_0)
  N(2) = 1. / 4. * (1. + c(0)) * (1. + c(1)); /// N3(q_0)
  N(3) = 1. / 4. * (1. - c(0)) * (1. + c(1)); /// N4(q_0)
}
/* -------------------------------------------------------------------------- */
template <>
template <class vector_type, class matrix_type>
inline void InterpolationElement<_itp_lagrange_quadrangle_4>::computeDNDS(
    const vector_type & c, matrix_type & dnds) {
  /**
   * @f[
   * dnds = \left(
   *          \begin{array}{cccc}
   *            \frac{\partial N1}{\partial \xi}  & \frac{\partial N2}{\partial
   * \xi}
   *               & \frac{\partial N3}{\partial \xi}  & \frac{\partial
   * N4}{\partial \xi}\\
   *            \frac{\partial N1}{\partial \eta} & \frac{\partial N2}{\partial
   * \eta}
   *               & \frac{\partial N3}{\partial \eta} & \frac{\partial
   * N4}{\partial \eta}
   *          \end{array}
   *        \right)
   * @f]
   */

  dnds(0, 0) = -1. / 4. * (1. - c(1));
  dnds(0, 1) = 1. / 4. * (1. - c(1));
  dnds(0, 2) = 1. / 4. * (1. + c(1));
  dnds(0, 3) = -1. / 4. * (1. + c(1));

  dnds(1, 0) = -1. / 4. * (1. - c(0));
  dnds(1, 1) = -1. / 4. * (1. + c(0));
  dnds(1, 2) = 1. / 4. * (1. + c(0));
  dnds(1, 3) = 1. / 4. * (1. - c(0));
}

/* -------------------------------------------------------------------------- */
template <>
inline void
InterpolationElement<_itp_lagrange_quadrangle_4>::computeSpecialJacobian(
    const Matrix<Real> & J, Real & jac) {
  Vector<Real> vprod(J.cols());
  Matrix<Real> Jt(J.transpose(), true);
  vprod.crossProduct(Jt(0), Jt(1));
  jac = vprod.norm();
}

/* -------------------------------------------------------------------------- */
template <>
inline Real
GeometricalElement<_gt_quadrangle_4>::getInradius(const Matrix<Real> & coord) {
  Vector<Real> u0 = coord(0);
  Vector<Real> u1 = coord(1);
  Vector<Real> u2 = coord(2);
  Vector<Real> u3 = coord(3);
  Real a = u0.distance(u1);
  Real b = u1.distance(u2);
  Real c = u2.distance(u3);
  Real d = u3.distance(u0);

  // Real septimetre = (a + b + c + d) / 2.;

  // Real p = Math::distance_2d(coord + 0, coord + 4);
  // Real q = Math::distance_2d(coord + 2, coord + 6);

  // Real area = sqrt(4*(p*p * q*q) - (a*a + b*b + c*c + d*d)*(a*a + c*c - b*b -
  // d*d)) / 4.;
  // Real h = sqrt(area);  // to get a length
  // Real h = area / septimetre;  // formula of inradius for circumscritable
  // quadrelateral
  Real h = std::min({a, b, c, d});

  return h;
}
} // namespace akantu
