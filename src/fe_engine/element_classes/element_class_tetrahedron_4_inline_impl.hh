/**
 * @file   element_class_tetrahedron_4_inline_impl.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jul 16 2010
 * @date last modification: Wed Oct 11 2017
 *
 * @brief  Specialization of the element_class class for the type _tetrahedron_4
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
     x (0,0,1,0)
     |`
     |  `  °             \zeta
     |    `    °       -
     |      `       x (0,0,0,1)
     |      q.`   -  '
     |         -`     '
     |     -       `   '
     | -             `  '
     x------------------x-----> \xi
 (1,0,0,0)           (0,1,0,0)
 @endverbatim
 *
 * @f{eqnarray*}{
 * N1 &=& 1 - \xi - \eta - \zeta \\
 * N2 &=& \xi \\
 * N3 &=& \eta \\
 * N4 &=& \zeta
 * @f}
 *
 * @f[
 * \xi_{q0} = 1/4 \qquad  \eta_{q0} = 1/4  \qquad  \zeta_{q0} = 1/4
 * @f]
 */
/* -------------------------------------------------------------------------- */
#include "element_class.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
/* -------------------------------------------------------------------------- */
AKANTU_DEFINE_ELEMENT_CLASS_PROPERTY(_tetrahedron_4, _gt_tetrahedron_4,
                                     _itp_lagrange_tetrahedron_4, _ek_regular,
                                     3, _git_tetrahedron, 1);

/* -------------------------------------------------------------------------- */
template <>
template <class vector_type>
inline void InterpolationElement<_itp_lagrange_tetrahedron_4>::computeShapes(
    const vector_type & natural_coords, vector_type & N) {

  Real c0 = 1 - natural_coords(0) - natural_coords(1) -
            natural_coords(2); /// @f$ c0 = 1 - \xi - \eta - \zeta @f$
  Real c1 = natural_coords(1); /// @f$ c1 = \xi @f$
  Real c2 = natural_coords(2); /// @f$ c2 = \eta @f$
  Real c3 = natural_coords(0); /// @f$ c3 = \zeta @f$

  N(0) = c0;
  N(1) = c1;
  N(2) = c2;
  N(3) = c3;
}
/* -------------------------------------------------------------------------- */
template <>
template <class vector_type, class matrix_type>
inline void InterpolationElement<_itp_lagrange_tetrahedron_4>::computeDNDS(
    __attribute__((unused)) const vector_type & natural_coords,
    matrix_type & dnds) {

  /**
   * @f[
   * dnds = \left(
   *          \begin{array}{cccccc}
   *            \frac{\partial N1}{\partial \xi}   & \frac{\partial N2}{\partial
   * \xi}
   *          & \frac{\partial N3}{\partial \xi}   & \frac{\partial N4}{\partial
   * \xi} \\
   *            \frac{\partial N1}{\partial \eta}  & \frac{\partial N2}{\partial
   * \eta}
   *          & \frac{\partial N3}{\partial \eta}  & \frac{\partial N4}{\partial
   * \eta} \\
   *            \frac{\partial N1}{\partial \zeta} & \frac{\partial N2}{\partial
   * \zeta}
   *          & \frac{\partial N3}{\partial \zeta} & \frac{\partial N4}{\partial
   * \zeta}
   *          \end{array}
   *        \right)
   * @f]
   */

  dnds(0, 0) = -1.;
  dnds(1, 0) = -1.;
  dnds(2, 0) = -1.;

  dnds(0, 1) = 0.;
  dnds(1, 1) = 1.;
  dnds(2, 1) = 0.;

  dnds(0, 2) = 0.;
  dnds(1, 2) = 0.;
  dnds(2, 2) = 1.;

  dnds(0, 3) = 1.;
  dnds(1, 3) = 0.;
  dnds(2, 3) = 0.;
}

/* -------------------------------------------------------------------------- */
template <>
inline Real
GeometricalElement<_gt_tetrahedron_4>::getInradius(const Matrix<Real> & coord) {
  return 2. * Math::tetrahedron_inradius(coord(0).storage(), coord(1).storage(),
                                         coord(2).storage(),
                                         coord(3).storage());
}
} // namespace akantu
