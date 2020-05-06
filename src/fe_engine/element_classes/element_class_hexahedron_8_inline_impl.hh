/**
 * @file   element_class_hexahedron_8_inline_impl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Peter Spijker <peter.spijker@epfl.ch>
 *
 * @date creation: Mon Mar 14 2011
 * @date last modification: Wed Oct 11 2017
 *
 * @brief  Specialization of the element_class class for the type _hexahedron_8
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
                   \zeta
                    ^
         (-1,1,1)   |     (1,1,1)
                7---|------6
               /|   |     /|
              / |   |    / |
   (-1,-1,1) 4----------5  | (1,-1,1)
             |  |   |   |  |
             |  |   |   |  |
             |  |   +---|-------> \xi
             |  |  /    |  |
   (-1,1,-1) |  3-/-----|--2 (1,1,-1)
             | / /      | /
             |/ /       |/
             0-/--------1
   (-1,-1,-1) /        (1,-1,-1)
             /
            \eta
 @endverbatim
 *
 * \f[
 * \begin{array}{llll}
 * N1 = (1 - \xi) (1 - \eta) (1 - \zeta) / 8
 *       & \frac{\partial N1}{\partial \xi}  = - (1 - \eta) (1 - \zeta) / 8
 *       & \frac{\partial N1}{\partial \eta} = - (1 - \xi) (1 - \zeta) / 8
 *       & \frac{\partial N1}{\partial \zeta} = - (1 - \xi) (1 - \eta) / 8 \\
 * N2 = (1 + \xi) (1 - \eta) (1 - \zeta) / 8
 *       & \frac{\partial N2}{\partial \xi}  = (1 - \eta) (1 - \zeta) / 8
 *       & \frac{\partial N2}{\partial \eta} = - (1 + \xi) (1 - \zeta) / 8
 *       & \frac{\partial N2}{\partial \zeta} = - (1 + \xi) (1 - \eta) / 8 \\
 * N3 = (1 + \xi) (1 + \eta) (1 - \zeta) / 8
 *       & \frac{\partial N3}{\partial \xi}  = (1 + \eta) (1 - \zeta) / 8
 *       & \frac{\partial N3}{\partial \eta} = (1 + \xi) (1 - \zeta) / 8
 *       & \frac{\partial N3}{\partial \zeta} = - (1 + \xi) (1 + \eta) / 8 \\
 * N4 = (1 - \xi) (1 + \eta) (1 - \zeta) / 8
 *       & \frac{\partial N4}{\partial \xi}  = - (1 + \eta) (1 - \zeta) / 8
 *       & \frac{\partial N4}{\partial \eta} = (1 - \xi) (1 - \zeta) / 8
 *       & \frac{\partial N4}{\partial \zeta} = - (1 - \xi) (1 + \eta) / 8 \\
 * N5 = (1 - \xi) (1 - \eta) (1 + \zeta) / 8
 *       & \frac{\partial N5}{\partial \xi}  = - (1 - \eta) (1 + \zeta) / 8
 *       & \frac{\partial N5}{\partial \eta} = - (1 - \xi) (1 + \zeta) / 8
 *       & \frac{\partial N5}{\partial \zeta} = (1 - \xi) (1 - \eta) / 8 \\
 * N6 = (1 + \xi) (1 - \eta) (1 + \zeta) / 8
 *       & \frac{\partial N6}{\partial \xi}  = (1 - \eta) (1 + \zeta) / 8
 *       & \frac{\partial N6}{\partial \eta} = - (1 + \xi) (1 + \zeta) / 8
 *       & \frac{\partial N6}{\partial \zeta} = (1 + \xi) (1 - \eta) / 8 \\
 * N7 = (1 + \xi) (1 + \eta) (1 + \zeta) / 8
 *       & \frac{\partial N7}{\partial \xi}  = (1 + \eta) (1 + \zeta) / 8
 *       & \frac{\partial N7}{\partial \eta} = (1 + \xi) (1 + \zeta) / 8
 *       & \frac{\partial N7}{\partial \zeta} = (1 + \xi) (1 + \eta) / 8 \\
 * N8 = (1 - \xi) (1 + \eta) (1 + \zeta) / 8
 *       & \frac{\partial N8}{\partial \xi}  = - (1 + \eta) (1 + \zeta) / 8
 *       & \frac{\partial N8}{\partial \eta} = (1 - \xi) (1 + \zeta) / 8
 *       & \frac{\partial N8}{\partial \zeta} = (1 - \xi) (1 + \eta) / 8 \\
 * \end{array}
 * \f]
 *
 * @f{eqnarray*}{
 * \xi_{q0}  &=& -1/\sqrt{3} \qquad  \eta_{q0} = -1/\sqrt{3} \qquad \zeta_{q0} =
 -1/\sqrt{3} \\
 * \xi_{q1}  &=&  1/\sqrt{3} \qquad  \eta_{q1} = -1/\sqrt{3} \qquad \zeta_{q1} =
 -1/\sqrt{3} \\
 * \xi_{q2}  &=&  1/\sqrt{3} \qquad  \eta_{q2} =  1/\sqrt{3} \qquad \zeta_{q2} =
 -1/\sqrt{3} \\
 * \xi_{q3}  &=& -1/\sqrt{3} \qquad  \eta_{q3} =  1/\sqrt{3} \qquad \zeta_{q3} =
 -1/\sqrt{3} \\
 * \xi_{q4}  &=& -1/\sqrt{3} \qquad  \eta_{q4} = -1/\sqrt{3} \qquad \zeta_{q4} =
 1/\sqrt{3} \\
 * \xi_{q5}  &=&  1/\sqrt{3} \qquad  \eta_{q5} = -1/\sqrt{3} \qquad \zeta_{q5} =
 1/\sqrt{3} \\
 * \xi_{q6}  &=&  1/\sqrt{3} \qquad  \eta_{q6} =  1/\sqrt{3} \qquad \zeta_{q6} =
 1/\sqrt{3} \\
 * \xi_{q7}  &=& -1/\sqrt{3} \qquad  \eta_{q7} =  1/\sqrt{3} \qquad \zeta_{q7} =
 1/\sqrt{3} \\
 * @f}
 */
/* -------------------------------------------------------------------------- */
#include "element_class.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
/* -------------------------------------------------------------------------- */
AKANTU_DEFINE_ELEMENT_CLASS_PROPERTY(_hexahedron_8, _gt_hexahedron_8,
                                     _itp_lagrange_hexahedron_8, _ek_regular, 3,
                                     _git_segment, 2);

/* -------------------------------------------------------------------------- */
template <>
template <class vector_type>
inline void InterpolationElement<_itp_lagrange_hexahedron_8>::computeShapes(
    const vector_type & c, vector_type & N) {
  /// Natural coordinates
  N(0) = .125 * (1 - c(0)) * (1 - c(1)) * (1 - c(2)); /// N1(q_0)
  N(1) = .125 * (1 + c(0)) * (1 - c(1)) * (1 - c(2)); /// N2(q_0)
  N(2) = .125 * (1 + c(0)) * (1 + c(1)) * (1 - c(2)); /// N3(q_0)
  N(3) = .125 * (1 - c(0)) * (1 + c(1)) * (1 - c(2)); /// N4(q_0)
  N(4) = .125 * (1 - c(0)) * (1 - c(1)) * (1 + c(2)); /// N5(q_0)
  N(5) = .125 * (1 + c(0)) * (1 - c(1)) * (1 + c(2)); /// N6(q_0)
  N(6) = .125 * (1 + c(0)) * (1 + c(1)) * (1 + c(2)); /// N7(q_0)
  N(7) = .125 * (1 - c(0)) * (1 + c(1)) * (1 + c(2)); /// N8(q_0)
}
/* -------------------------------------------------------------------------- */
template <>
template <class vector_type, class matrix_type>
inline void InterpolationElement<_itp_lagrange_hexahedron_8>::computeDNDS(
    const vector_type & c, matrix_type & dnds) {
  /**
   * @f[
   * dnds = \left(
   *          \begin{array}{cccccccc}
   *            \frac{\partial N1}{\partial \xi}  & \frac{\partial N2}{\partial
   * \xi}
   *               & \frac{\partial N3}{\partial \xi}  & \frac{\partial
   * N4}{\partial \xi}
   *               & \frac{\partial N5}{\partial \xi}  & \frac{\partial
   * N6}{\partial \xi}
   *               & \frac{\partial N7}{\partial \xi}  & \frac{\partial
   * N8}{\partial \xi}\\
   *            \frac{\partial N1}{\partial \eta} & \frac{\partial N2}{\partial
   * \eta}
   *               & \frac{\partial N3}{\partial \eta} & \frac{\partial
   * N4}{\partial \eta}
   *               & \frac{\partial N5}{\partial \eta} & \frac{\partial
   * N6}{\partial \eta}
   *               & \frac{\partial N7}{\partial \eta} & \frac{\partial
   * N8}{\partial \eta}\\
   *            \frac{\partial N1}{\partial \zeta} & \frac{\partial N2}{\partial
   * \zeta}
   *               & \frac{\partial N3}{\partial \zeta} & \frac{\partial
   * N4}{\partial \zeta}
   *               & \frac{\partial N5}{\partial \zeta} & \frac{\partial
   * N6}{\partial \zeta}
   *               & \frac{\partial N7}{\partial \zeta} & \frac{\partial
   * N8}{\partial \zeta}
   *          \end{array}
   *        \right)
   * @f]
   */
  dnds(0, 0) = -.125 * (1 - c(1)) * (1 - c(2));
  dnds(0, 1) = .125 * (1 - c(1)) * (1 - c(2));
  dnds(0, 2) = .125 * (1 + c(1)) * (1 - c(2));
  dnds(0, 3) = -.125 * (1 + c(1)) * (1 - c(2));
  dnds(0, 4) = -.125 * (1 - c(1)) * (1 + c(2));
  ;
  dnds(0, 5) = .125 * (1 - c(1)) * (1 + c(2));
  ;
  dnds(0, 6) = .125 * (1 + c(1)) * (1 + c(2));
  ;
  dnds(0, 7) = -.125 * (1 + c(1)) * (1 + c(2));
  ;

  dnds(1, 0) = -.125 * (1 - c(0)) * (1 - c(2));
  ;
  dnds(1, 1) = -.125 * (1 + c(0)) * (1 - c(2));
  ;
  dnds(1, 2) = .125 * (1 + c(0)) * (1 - c(2));
  ;
  dnds(1, 3) = .125 * (1 - c(0)) * (1 - c(2));
  ;
  dnds(1, 4) = -.125 * (1 - c(0)) * (1 + c(2));
  ;
  dnds(1, 5) = -.125 * (1 + c(0)) * (1 + c(2));
  ;
  dnds(1, 6) = .125 * (1 + c(0)) * (1 + c(2));
  ;
  dnds(1, 7) = .125 * (1 - c(0)) * (1 + c(2));
  ;

  dnds(2, 0) = -.125 * (1 - c(0)) * (1 - c(1));
  ;
  dnds(2, 1) = -.125 * (1 + c(0)) * (1 - c(1));
  ;
  dnds(2, 2) = -.125 * (1 + c(0)) * (1 + c(1));
  ;
  dnds(2, 3) = -.125 * (1 - c(0)) * (1 + c(1));
  ;
  dnds(2, 4) = .125 * (1 - c(0)) * (1 - c(1));
  ;
  dnds(2, 5) = .125 * (1 + c(0)) * (1 - c(1));
  ;
  dnds(2, 6) = .125 * (1 + c(0)) * (1 + c(1));
  ;
  dnds(2, 7) = .125 * (1 - c(0)) * (1 + c(1));
  ;
}

/* -------------------------------------------------------------------------- */
template <>
inline Real
GeometricalElement<_gt_hexahedron_8>::getInradius(const Matrix<Real> & coord) {
  Vector<Real> u0 = coord(0);
  Vector<Real> u1 = coord(1);
  Vector<Real> u2 = coord(2);
  Vector<Real> u3 = coord(3);
  Vector<Real> u4 = coord(4);
  Vector<Real> u5 = coord(5);
  Vector<Real> u6 = coord(6);
  Vector<Real> u7 = coord(7);

  Real a = u0.distance(u1);
  Real b = u1.distance(u2);
  Real c = u2.distance(u3);
  Real d = u3.distance(u0);
  Real e = u0.distance(u4);
  Real f = u1.distance(u5);
  Real g = u2.distance(u6);
  Real h = u3.distance(u7);
  Real i = u4.distance(u5);
  Real j = u5.distance(u6);
  Real k = u6.distance(u7);
  Real l = u7.distance(u4);

  Real p = std::min({a, b, c, d, e, f, g, h, i, j, k, l});

  return p;
}
} // namespace akantu
