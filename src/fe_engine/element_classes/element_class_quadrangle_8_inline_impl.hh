/**
 * @file   element_class_quadrangle_8_inline_impl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed May 18 2011
 * @date last modification: Wed Oct 11 2017
 *
 * @brief  Specialization of the ElementClass for the _quadrangle_8
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
   (-1,1)    (0,1)   (1,1)
       x-------x-------x
       |       |       |
       |       |       |
       |       |       |
 (-1,0)|       |       |(1,0)
   ----x---------------X----->  \xi
       |       |       |
       |       |       |
       |       |       |
       |       |       |
       x-------x-------x
   (-1,-1)   (0,-1)  (1,-1)
               |
 @endverbatim
 *
 * @f[
 * \begin{array}{lll}
 * N1 = (1 - \xi) (1 - \eta)(- 1 - \xi - \eta) / 4
 *       & \frac{\partial N1}{\partial \xi}  = (1 - \eta)(2 \xi + \eta) / 4
 *       & \frac{\partial N1}{\partial \eta} = (1 - \xi)(\xi + 2 \eta) / 4 \\
 * N2 = (1 + \xi) (1 - \eta)(- 1 + \xi - \eta) / 4 \\
 *       & \frac{\partial N2}{\partial \xi}  = (1 - \eta)(2 \xi - \eta) / 4
 *       & \frac{\partial N2}{\partial \eta} = - (1 + \xi)(\xi - 2 \eta) / 4 \\
 * N3 = (1 + \xi) (1 + \eta)(- 1 + \xi + \eta) / 4 \\
 *       & \frac{\partial N3}{\partial \xi}  = (1 + \eta)(2 \xi + \eta) / 4
 *       & \frac{\partial N3}{\partial \eta} = (1 + \xi)(\xi + 2 \eta) / 4 \\
 * N4 = (1 - \xi) (1 + \eta)(- 1 - \xi + \eta) / 4
 *       & \frac{\partial N4}{\partial \xi}  = (1 + \eta)(2 \xi - \eta) / 4
 *       & \frac{\partial N4}{\partial \eta} = - (1 - \xi)(\xi - 2 \eta) / 4 \\
 * N5 = (1 - \xi^2) (1 - \eta) / 2
 *       & \frac{\partial N1}{\partial \xi}  = - \xi (1 - \eta)
 *       & \frac{\partial N1}{\partial \eta} = - (1 - \xi^2) / 2  \\
 * N6 = (1 + \xi) (1 - \eta^2) / 2 \\
 *       & \frac{\partial N2}{\partial \xi}  = (1 - \eta^2) / 2
 *       & \frac{\partial N2}{\partial \eta} = - \eta (1 + \xi) \\
 * N7 = (1 - \xi^2) (1 + \eta) / 2 \\
 *       & \frac{\partial N3}{\partial \xi}  = - \xi (1 + \eta)
 *       & \frac{\partial N3}{\partial \eta} = (1 - \xi^2) / 2 \\
 * N8 = (1 - \xi) (1 - \eta^2) / 2
 *       & \frac{\partial N4}{\partial \xi}  = - (1 - \eta^2) / 2
 *       & \frac{\partial N4}{\partial \eta} = - \eta (1 - \xi) \\
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
AKANTU_DEFINE_ELEMENT_CLASS_PROPERTY(_quadrangle_8, _gt_quadrangle_8,
                                     _itp_serendip_quadrangle_8, _ek_regular, 2,
                                     _git_segment, 3);

/* -------------------------------------------------------------------------- */
template <>
template <class vector_type>
inline void InterpolationElement<_itp_serendip_quadrangle_8>::computeShapes(
    const vector_type & c, vector_type & N) {

  /// Natural coordinates
  const Real xi = c(0);
  const Real eta = c(1);

  N(0) = .25 * (1 - xi) * (1 - eta) * (-1 - xi - eta);
  N(1) = .25 * (1 + xi) * (1 - eta) * (-1 + xi - eta);
  N(2) = .25 * (1 + xi) * (1 + eta) * (-1 + xi + eta);
  N(3) = .25 * (1 - xi) * (1 + eta) * (-1 - xi + eta);
  N(4) = .5 * (1 - xi * xi) * (1 - eta);
  N(5) = .5 * (1 + xi) * (1 - eta * eta);
  N(6) = .5 * (1 - xi * xi) * (1 + eta);
  N(7) = .5 * (1 - xi) * (1 - eta * eta);
}

/* -------------------------------------------------------------------------- */
template <>
template <class vector_type, class matrix_type>
inline void InterpolationElement<_itp_serendip_quadrangle_8>::computeDNDS(
    const vector_type & c, matrix_type & dnds) {

  const Real xi = c(0);
  const Real eta = c(1);

  /// dN/dxi
  dnds(0, 0) = .25 * (1 - eta) * (2 * xi + eta);
  dnds(0, 1) = .25 * (1 - eta) * (2 * xi - eta);
  dnds(0, 2) = .25 * (1 + eta) * (2 * xi + eta);
  dnds(0, 3) = .25 * (1 + eta) * (2 * xi - eta);
  dnds(0, 4) = -xi * (1 - eta);
  dnds(0, 5) = .5 * (1 - eta * eta);
  dnds(0, 6) = -xi * (1 + eta);
  dnds(0, 7) = -.5 * (1 - eta * eta);

  /// dN/deta
  dnds(1, 0) = .25 * (1 - xi) * (2 * eta + xi);
  dnds(1, 1) = .25 * (1 + xi) * (2 * eta - xi);
  dnds(1, 2) = .25 * (1 + xi) * (2 * eta + xi);
  dnds(1, 3) = .25 * (1 - xi) * (2 * eta - xi);
  dnds(1, 4) = -.5 * (1 - xi * xi);
  dnds(1, 5) = -eta * (1 + xi);
  dnds(1, 6) = .5 * (1 - xi * xi);
  dnds(1, 7) = -eta * (1 - xi);
}

/* -------------------------------------------------------------------------- */
template <>
inline Real
GeometricalElement<_gt_quadrangle_8>::getInradius(const Matrix<Real> & coord) {
  Vector<Real> u0 = coord(0);
  Vector<Real> u1 = coord(1);
  Vector<Real> u2 = coord(2);
  Vector<Real> u3 = coord(3);
  Vector<Real> u4 = coord(4);
  Vector<Real> u5 = coord(5);
  Vector<Real> u6 = coord(6);
  Vector<Real> u7 = coord(7);

  auto a = u0.distance(u4);
  auto b = u4.distance(u1);
  auto h = std::min(a, b);

  a = u1.distance(u5);
  b = u5.distance(u2);
  h = std::min(h, std::min(a, b));

  a = u2.distance(u6);
  b = u6.distance(u3);
  h = std::min(h, std::min(a, b));

  a = u3.distance(u7);
  b = u7.distance(u0);
  h = std::min(h, std::min(a, b));

  return h;
}

/* -------------------------------------------------------------------------- */
template <>
inline void
InterpolationElement<_itp_serendip_quadrangle_8>::computeSpecialJacobian(
    const Matrix<Real> & J, Real & jac) {
  Vector<Real> vprod(J.cols());
  Matrix<Real> Jt(J.transpose(), true);
  vprod.crossProduct(Jt(0), Jt(1));
  jac = vprod.norm();
}
} // namespace akantu
