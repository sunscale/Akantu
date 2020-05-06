/**
 * @file   element_class_igfem_triangle_5_inline_impl.hh
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Tue May  5 14:57:51 2015
 *
 * @brief Specialization of the element_class class for the type
 * _igfem_triangle_5
 *
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * @section DESCRIPTION
 *
 * @verbatim
       \eta
     ^
     |
     x (0,0,1)
     |`
     |  `
     |` q `
     |  ` ° `
     x--------x----->  \xi
    (1,0,0)      (0,1,0)
 @endverbatim
 *
 * @subsection shapes shape functions
 * Parent:
 * @f{eqnarray*}{
 * N1 &=& 1 - \xi - \eta \\
 * N2 &=& \xi \\
 * N3 &=& \eta
 * @f}
 * Sub 1:
 * @f{eqnarray*}{
 * N2 &=& \xi \\
 * N3 &=& \eta
 * @f}
 * Sub 2:
 * @f[
 * \begin{array}{lll}
 * N3 = (1 + \xi) (1 + \eta) / 4 \\
 *       & \frac{\partial N3}{\partial \xi}  = (1 + \eta) / 4
 *       & \frac{\partial N3}{\partial \eta} = (1 + \xi) / 4 \\
 * N4 = (1 - \xi) (1 + \eta) / 4
 *       & \frac{\partial N4}{\partial \xi}  = - (1 + \eta) / 4
 *       & \frac{\partial N4}{\partial \eta} = (1 - \xi) / 4 \\
 * \end{array}
 * @f]
 *
 * @subsection quad_points Position of quadrature points
 * @f{eqnarray*}{
 * \xi_{q0}  &=& 1/3 \qquad  \eta_{q0} = 1/3
 * @f}
 */

/* -------------------------------------------------------------------------- */
#include "element_class_igfem.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
AKANTU_DEFINE_IGFEM_ELEMENT_CLASS_PROPERTY(_igfem_triangle_5,
                                           _gt_igfem_triangle_5,
                                           _itp_igfem_triangle_5, _triangle_3,
                                           _triangle_3, _quadrangle_4,
                                           _ek_igfem, 2, 1);

/* -------------------------------------------------------------------------- */
template <>
inline UInt ElementClass<_igfem_triangle_5>::getOrientation(
    const Vector<bool> & is_inside) {
  UInt sub_el_is_inside = 0;
  if (is_inside(0)) {
    sub_el_is_inside = 0;

    AKANTU_DEBUG_ASSERT(!(is_inside(1) || is_inside(2)),
                        "orientation not determinable");
  } else {
    sub_el_is_inside = 1;
    AKANTU_DEBUG_ASSERT((is_inside(2) && is_inside(2)),
                        "orientation not determinable");
  }
  return sub_el_is_inside;
}
} // namespace akantu
