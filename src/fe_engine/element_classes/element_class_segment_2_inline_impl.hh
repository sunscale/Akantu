/**
 * @file   element_class_segment_2_inline_impl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jul 16 2010
 * @date last modification: Wed Oct 11 2017
 *
 * @brief  Specialization of the element_class class for the type _segment_2
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
              q
   --x--------|--------x---> x
    -1        0        1
 @endverbatim
 *
 * @f{eqnarray*}{
 * w_1(x) &=& 1/2(1 - x) \\
 * w_2(x) &=& 1/2(1 + x)
 * @f}
 *
 * @f{eqnarray*}{
 * x_{q}  &=& 0
 * @f}
 */
/* -------------------------------------------------------------------------- */
#include "element_class.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
/* -------------------------------------------------------------------------- */
AKANTU_DEFINE_ELEMENT_CLASS_PROPERTY(_segment_2, _gt_segment_2,
                                     _itp_lagrange_segment_2, _ek_regular, 1,
                                     _git_segment, 1);

/* -------------------------------------------------------------------------- */
template <>
template <class vector_type>
inline void InterpolationElement<_itp_lagrange_segment_2>::computeShapes(
    const vector_type & natural_coords, vector_type & N) {

  /// natural coordinate
  Real c = natural_coords(0);
  /// shape functions
  N(0) = 0.5 * (1 - c);
  N(1) = 0.5 * (1 + c);
}
/* -------------------------------------------------------------------------- */

template <>
template <class vector_type, class matrix_type>
inline void InterpolationElement<_itp_lagrange_segment_2>::computeDNDS(
    __attribute__((unused)) const vector_type & natural_coords,
    matrix_type & dnds) {

  /// dN1/de
  dnds(0, 0) = -.5;
  /// dN2/de
  dnds(0, 1) = .5;
}

/* -------------------------------------------------------------------------- */
template <>
inline void
InterpolationElement<_itp_lagrange_segment_2>::computeSpecialJacobian(
    const Matrix<Real> & dxds, Real & jac) {
  jac = dxds.norm<L_2>();
}

/* -------------------------------------------------------------------------- */
template <>
inline Real
GeometricalElement<_gt_segment_2>::getInradius(const Matrix<Real> & coord) {
  return std::abs(coord(0, 0) - coord(0, 1));
}

// /* --------------------------------------------------------------------------
// */
// template<> inline bool ElementClass<_segment_2>::contains(const Vector<Real>
// & natural_coords) {
//   if (natural_coords(0) < -1.) return false;
//   if (natural_coords(0) > 1.) return false;
//   return true;
// }

/* -------------------------------------------------------------------------- */
} // namespace akantu
