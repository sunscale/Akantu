/**
 * @file   element_class_segment_3_inline_impl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jul 16 2010
 * @date last modification: Wed Oct 11 2017
 *
 * @brief  Specialization of the element_class class for the type _segment_3
 *
 * \section LICENSE
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
         -1         0         1
     -----x---------x---------x-----> x
          1         3         2
 @endverbatim
 *
 *
 * @f[
 * \begin{array}{lll}
 *   x_{1}  = -1   &  x_{2} = 1 & x_{3} = 0
 * \end{array}
 * @f]
 *
 * @f[
 * \begin{array}{ll}
 *   w_1(x) = \frac{x}{2}(x - 1) & w'_1(x) = x - \frac{1}{2}\\
 *   w_2(x) = \frac{x}{2}(x + 1) & w'_2(x) = x + \frac{1}{2}\\
 *   w_3(x) =  1-x^2 & w'_3(x) = -2x
 * \end{array}
 * @f]
 *
 * @f[
 * \begin{array}{ll}
 * x_{q1}  = -1/\sqrt{3} & x_{q2} = 1/\sqrt{3}
 * \end{array}
 * @f]
 */
/* -------------------------------------------------------------------------- */
#include "element_class.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
/* -------------------------------------------------------------------------- */
AKANTU_DEFINE_ELEMENT_CLASS_PROPERTY(_segment_3, _gt_segment_3,
                                     _itp_lagrange_segment_3, _ek_regular, 1,
                                     _git_segment, 2);

/* -------------------------------------------------------------------------- */
template <>
template <class vector_type>
inline void InterpolationElement<_itp_lagrange_segment_3>::computeShapes(
    const vector_type & natural_coords, vector_type & N) {

  Real c = natural_coords(0);
  N(0) = (c - 1) * c / 2;
  N(1) = (c + 1) * c / 2;
  N(2) = 1 - c * c;
}
/* -------------------------------------------------------------------------- */
template <>
template <class vector_type, class matrix_type>
inline void InterpolationElement<_itp_lagrange_segment_3>::computeDNDS(
    const vector_type & natural_coords, matrix_type & dnds) {

  Real c = natural_coords(0);
  dnds(0, 0) = c - .5;
  dnds(0, 1) = c + .5;
  dnds(0, 2) = -2 * c;
}

/* -------------------------------------------------------------------------- */
template <>
inline void
InterpolationElement<_itp_lagrange_segment_3>::computeSpecialJacobian(
    const Matrix<Real> & dxds, Real & jac) {
  jac = Math::norm2(dxds.storage());
}

/* -------------------------------------------------------------------------- */
template <>
inline Real
GeometricalElement<_gt_segment_3>::getInradius(const Matrix<Real> & coord) {
  Real dist1 = std::abs(coord(0, 0) - coord(0, 1));
  Real dist2 = std::abs(coord(0, 1) - coord(0, 2));
  return std::min(dist1, dist2);
}
} // namespace akantu
