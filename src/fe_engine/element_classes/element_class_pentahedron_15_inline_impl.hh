/**
 * @file   element_class_pentahedron_15_inline_impl.hh
 *
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 * @author Sacha Laffely <sacha.laffely@epfl.ch>
 * @author Damien Scantamburlo <damien.scantamburlo@epfl.ch>
 *
 * @date creation: Tue Mar 31 2015
 * @date last modification: Thu Dec 28 2017
 *
 * @brief  Specialization of the element_class class for the type
 * _pentahedron_15
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
 * \verbatim
             z
             ^
             |
             |
             |  1
             | /|\
             |/ | \
             10 7  6
            /   |   \
           /    |    \
           4    2--8--0
          | \  /      /
          |  \11     /
          13  12    9---------->y
          | /  \   /
          |/    \ /
          5--14--3
         /
        /
       /
      v
     x
\endverbatim
       x   y    z
* N0  -1   1    0
* N1  -1   0    1
* N2  -1   0    0
* N3   1   1    0
* N4   1   0    1
* N5   1   0    0
* N6  -1   0.5  0.5
* N7  -1   0    0.5
* N8  -1   0.5  0
* N9   0   1    0
* N10  0   0    1
* N11  0   0    0
* N12  1   0.5  0.5
* N13  1   0    0.5
* N14  1   0.5  0
 */

/* -------------------------------------------------------------------------- */
#include "element_class.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
/* -------------------------------------------------------------------------- */
AKANTU_DEFINE_ELEMENT_CLASS_PROPERTY(_pentahedron_15, _gt_pentahedron_15,
                                     _itp_lagrange_pentahedron_15, _ek_regular,
                                     3, _git_pentahedron, 2);

/* -------------------------------------------------------------------------- */
template <>
template <class vector_type>
inline void InterpolationElement<_itp_lagrange_pentahedron_15>::computeShapes(
    const vector_type & c, vector_type & N) {
  auto & x = c(0);
  auto & y = c(1);
  auto & z = c(2);

  // Shape Functions, Natural coordinates
  N(0) = 0.5 * y * (1 - x) * (2 * y - 2 - x);
  N(1) = 0.5 * z * (1 - x) * (2 * z - 2 - x);
  N(2) = 0.5 * (x - 1) * (1 - y - z) * (x + 2 * y + 2 * z);
  N(3) = 0.5 * y * (1 + x) * (2 * y - 2 + x);
  N(4) = 0.5 * z * (1 + x) * (2 * z - 2 + x);
  N(5) = 0.5 * (-x - 1) * (1 - y - z) * (-x + 2 * y + 2 * z);
  N(6) = 2.0 * y * z * (1 - x);
  N(7) = 2.0 * z * (1 - y - z) * (1 - x);
  N(8) = 2.0 * y * (1 - x) * (1 - y - z);
  N(9) = y * (1 - x * x);
  N(10) = z * (1 - x * x);
  N(11) = (1 - y - z) * (1 - x * x);
  N(12) = 2.0 * y * z * (1 + x);
  N(13) = 2.0 * z * (1 - y - z) * (1 + x);
  N(14) = 2.0 * y * (1 - y - z) * (1 + x);
}

/* -------------------------------------------------------------------------- */
template <>
template <class vector_type, class matrix_type>
inline void InterpolationElement<_itp_lagrange_pentahedron_15>::computeDNDS(
    const vector_type & c, matrix_type & dnds) {
  auto & x = c(0);
  auto & y = c(1);
  auto & z = c(2);

  // ddx
  dnds(0, 0) = 0.5 * y * (2 * x - 2 * y + 1);
  dnds(0, 1) = 0.5 * z * (2 * x - 2 * z + 1);
  dnds(0, 2) = -0.5 * (2 * x + 2 * y + 2 * z - 1) * (y + z - 1);
  dnds(0, 3) = 0.5 * y * (2 * x + 2 * y - 1);
  dnds(0, 4) = 0.5 * z * (2 * x + 2 * z - 1);
  dnds(0, 5) = -0.5 * (y + z - 1) * (2 * x - 2 * y - 2 * z + 1);
  dnds(0, 6) = -2.0 * y * z;
  dnds(0, 7) = 2.0 * z * (y + z - 1);
  dnds(0, 8) = 2.0 * y * (y + z - 1);
  dnds(0, 9) = -2.0 * x * y;
  dnds(0, 10) = -2.0 * x * z;
  dnds(0, 11) = 2.0 * x * (y + z - 1);
  dnds(0, 12) = 2.0 * y * z;
  dnds(0, 13) = -2.0 * z * (y + z - 1);
  dnds(0, 14) = -2.0 * y * (y + z - 1);

  // ddy
  dnds(1, 0) = -0.5 * (x - 1) * (4 * y - x - 2);
  dnds(1, 1) = 0.0;
  dnds(1, 2) = -0.5 * (x - 1) * (4 * y + x + 2 * (2 * z - 1));
  dnds(1, 3) = 0.5 * (x + 1) * (4 * y + x - 2);
  dnds(1, 4) = 0.0;
  dnds(1, 5) = 0.5 * (x + 1) * (4 * y - x + 2 * (2 * z - 1));
  dnds(1, 6) = -2.0 * (x - 1) * z;
  dnds(1, 7) = 2.0 * z * (x - 1);
  dnds(1, 8) = 2.0 * (2 * y + z - 1) * (x - 1);
  dnds(1, 9) = -(x * x - 1);
  dnds(1, 10) = 0.0;
  dnds(1, 11) = (x * x - 1);
  dnds(1, 12) = 2.0 * z * (x + 1);
  dnds(1, 13) = -2.0 * z * (x + 1);
  dnds(1, 14) = -2.0 * (2 * y + z - 1) * (x + 1);

  // ddz
  dnds(2, 0) = 0.0;
  dnds(2, 1) = -0.5 * (x - 1) * (4 * z - x - 2);
  dnds(2, 2) = -0.5 * (x - 1) * (4 * z + x + 2 * (2 * y - 1));
  dnds(2, 3) = 0.0;
  dnds(2, 4) = 0.5 * (x + 1) * (4 * z + x - 2);
  dnds(2, 5) = 0.5 * (x + 1) * (4 * z - x + 2 * (2 * y - 1));
  dnds(2, 6) = -2.0 * (x - 1) * y;
  dnds(2, 7) = 2.0 * (x - 1) * (2 * z + y - 1);
  dnds(2, 8) = 2.0 * y * (x - 1);
  dnds(2, 9) = 0.0;
  dnds(2, 10) = -(x * x - 1);
  dnds(2, 11) = (x * x - 1);
  dnds(2, 12) = 2.0 * (x + 1) * y;
  dnds(2, 13) = -2.0 * (x + 1) * (2 * z + y - 1);
  dnds(2, 14) = -2.0 * (x + 1) * y;
}

/* -------------------------------------------------------------------------- */
template <>
inline Real GeometricalElement<_gt_pentahedron_15>::getInradius(
    const Matrix<Real> & coord) {
  return GeometricalElement<_gt_pentahedron_6>::getInradius(coord) * 0.5;
}
} // namespace akantu
