/**
 * @file   element_class_pentahedron_6_inline_impl.hh
 *
 * @author Marion Estelle Chambart <mchambart@stucky.ch>
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 * @author Thomas Menouillard <tmenouillard@stucky.ch>
 *
 * @date creation: Mon Mar 14 2011
 * @date last modification: Wed Oct 11 2017
 *
 * @brief  Specialization of the element_class class for the type _pentahedron_6
 *
 *
 * Copyright (©) 2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * @verbatim
             /z
             |
             |
             |  1
             | /|\
             |/ | \
             /  |  \
            /   |   \
           /    |    \
           4    2-----0
          | \  /      /
          |  \/      /
          |   \     /----------/y
          | /  \   /
          |/    \ /
          5---.--3
         /
        /
       /
      \x
       x   y    z
* N0  -1   1    0
* N1  -1   0    1
* N2  -1   0    0
* N3   1   1    0
* N4   1   0    1
* N5   1   0    0
 \endverbatim
 */
/* -------------------------------------------------------------------------- */
#include "element_class.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
/* -------------------------------------------------------------------------- */
AKANTU_DEFINE_ELEMENT_CLASS_PROPERTY(_pentahedron_6, _gt_pentahedron_6,
                                     _itp_lagrange_pentahedron_6, _ek_regular,
                                     3, _git_pentahedron, 1);

/* -------------------------------------------------------------------------- */
template <>
template <class vector_type>
inline void InterpolationElement<_itp_lagrange_pentahedron_6>::computeShapes(
    const vector_type & c, vector_type & N) {
  /// Natural coordinates
  N(0) = 0.5 * c(1) * (1 - c(0));              // N1(q)
  N(1) = 0.5 * c(2) * (1 - c(0));              // N2(q)
  N(2) = 0.5 * (1 - c(1) - c(2)) * (1 - c(0)); // N3(q)
  N(3) = 0.5 * c(1) * (1 + c(0));              // N4(q)
  N(4) = 0.5 * c(2) * (1 + c(0));              // N5(q)
  N(5) = 0.5 * (1 - c(1) - c(2)) * (1 + c(0)); // N6(q)
}
/* -------------------------------------------------------------------------- */
template <>
template <class vector_type, class matrix_type>
inline void InterpolationElement<_itp_lagrange_pentahedron_6>::computeDNDS(
    const vector_type & c, matrix_type & dnds) {
  dnds(0, 0) = -0.5 * c(1);
  dnds(0, 1) = -0.5 * c(2);
  dnds(0, 2) = -0.5 * (1 - c(1) - c(2));
  dnds(0, 3) =  0.5 * c(1);
  dnds(0, 4) =  0.5 * c(2);
  dnds(0, 5) =  0.5 * (1 - c(1) - c(2));

  dnds(1, 0) = 0.5 * (1 - c(0));
  dnds(1, 1) = 0.0;
  dnds(1, 2) = -0.5 * (1 - c(0));
  dnds(1, 3) = 0.5 * (1 + c(0));
  dnds(1, 4) = 0.0;
  dnds(1, 5) = -0.5 * (1 + c(0));

  dnds(2, 0) = 0.0;
  dnds(2, 1) = 0.5 * (1 - c(0));
  dnds(2, 2) = -0.5 * (1 - c(0));
  dnds(2, 3) = 0.0;
  dnds(2, 4) = 0.5 * (1 + c(0));
  dnds(2, 5) = -0.5 * (1 + c(0));
}

/* -------------------------------------------------------------------------- */
// I have to duplicate this code since the Real * coords do not know their size
// in the Math module.
// If later we use eigen or Vector to implement this function
// there should be only one function in akantu::Math
// -> this is temporary for the release deadline which was so extended

inline Real triangle_inradius(const Real * coord1, const Real * coord2,
                              const Real * coord3) {
  /**
   * @f{eqnarray*}{
   * r &=& A / s \\
   * A &=& 1/4 * \sqrt{(a + b + c) * (a - b + c) * (a + b - c) (-a + b + c)} \\
   * s &=& \frac{a + b + c}{2}
   * @f}
   */

  auto a = Math::distance_3d(coord1, coord2);
  auto b = Math::distance_3d(coord2, coord3);
  auto c = Math::distance_3d(coord1, coord3);

  auto s = (a + b + c) * 0.5;

  return std::sqrt((s - a) * (s - b) * (s - c) / s);
}
/* -------------------------------------------------------------------------- */
template <>
inline Real
GeometricalElement<_gt_pentahedron_6>::getInradius(const Matrix<Real> & coord) {
  Vector<Real> u0 = coord(0);
  Vector<Real> u1 = coord(1);
  Vector<Real> u2 = coord(2);
  Vector<Real> u3 = coord(3);
  Vector<Real> u4 = coord(4);
  Vector<Real> u5 = coord(5);

  auto inradius_triangle_1 =
      triangle_inradius(u0.storage(), u1.storage(), u2.storage());

  auto inradius_triangle_2 =
      triangle_inradius(u3.storage(), u4.storage(), u5.storage());

  auto d1 = u3.distance(u0) * 0.5;
  auto d2 = u5.distance(u2) * 0.5;
  auto d3 = u4.distance(u1) * 0.5;
  auto p =
      2. * std::min({inradius_triangle_1, inradius_triangle_2, d1, d2, d3});

  return p;
}
} // namespace akantu
