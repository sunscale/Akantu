/**
 * @file   element_class_hexahedron_20_inline_impl.cc
 *
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 * @author Sacha Laffely <sacha.laffely@epfl.ch>
 * @author Damien Scantamburlo <damien.scantamburlo@epfl.ch>
 *
 * @date creation: Tue Mar 31 2015
 * @date last modification: Thu Jul 16 2015
 *
 * @brief  Specialization of the element_class class for the type _hexahedron_20
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
 * @section DESCRIPTION
 *
 * @verbatim
                                   \y
                         \z       /
                         |       /
                   7-----|18--------6
                  /|     |     /   /|
                 / |     |    /   / |
               19  |     |   /  17  |
               /  15     |  /   /   14
              /    |     | /   /    |
             4-------16---/---5     |
             |     |     +----|------------\x
             |     3-------10-|-----2
             |    /           |    /
            12   /           13   /
             |  11            |  9
             | /              | /
             |/               |/
             0--------8-------1
       x   y    z
* N0  -1  -1   -1
* N1   1  -1   -1
* N2   1   1   -1
* N3  -1   1   -1
* N4  -1  -1    1
* N5   1  -1    1
* N6   1   1    1
* N7  -1   1    1
* N8   0  -1   -1
* N9   1   0   -1
* N10  0   1   -1
* N11 -1   0   -1
* N12 -1  -1    0
* N13  1  -1    0
* N14  1   1    0
* N15 -1   1    0
* N16  0  -1    1
* N17  1   0    1
* N18  0   1    1
* N19 -1   0    1
 */

/* -------------------------------------------------------------------------- */
AKANTU_DEFINE_ELEMENT_CLASS_PROPERTY(_hexahedron_20, _gt_hexahedron_20,
                                     _itp_serendip_hexahedron_20, _ek_regular,
                                     3, _git_segment, 3);

/* -------------------------------------------------------------------------- */
template <>
template <class vector_type>
inline void InterpolationElement<_itp_serendip_hexahedron_20>::computeShapes(
    const vector_type & c, vector_type & N) {

  // Shape function , Natural coordinates
  N(0) =
      0.125 * (1 - c(0)) * (1 - c(1)) * (1 - c(2)) * (-2 - c(0) - c(1) - c(2));
  N(1) =
      0.125 * (1 + c(0)) * (1 - c(1)) * (1 - c(2)) * (-2 + c(0) - c(1) - c(2));
  N(2) =
      0.125 * (1 + c(0)) * (1 + c(1)) * (1 - c(2)) * (-2 + c(0) + c(1) - c(2));
  N(3) =
      0.125 * (1 - c(0)) * (1 + c(1)) * (1 - c(2)) * (-2 - c(0) + c(1) - c(2));
  N(4) =
      0.125 * (1 - c(0)) * (1 - c(1)) * (1 + c(2)) * (-2 - c(0) - c(1) + c(2));
  N(5) =
      0.125 * (1 + c(0)) * (1 - c(1)) * (1 + c(2)) * (-2 + c(0) - c(1) + c(2));
  N(6) =
      0.125 * (1 + c(0)) * (1 + c(1)) * (1 + c(2)) * (-2 + c(0) + c(1) + c(2));
  N(7) =
      0.125 * (1 - c(0)) * (1 + c(1)) * (1 + c(2)) * (-2 - c(0) + c(1) + c(2));
  N(8) = 0.25 * (1 - c(0) * c(0)) * (1 - c(1)) * (1 - c(2));
  N(9) = 0.25 * (1 - c(1) * c(1)) * (1 + c(0)) * (1 - c(2));
  N(10) = 0.25 * (1 - c(0) * c(0)) * (1 + c(1)) * (1 - c(2));
  N(11) = 0.25 * (1 - c(1) * c(1)) * (1 - c(0)) * (1 - c(2));
  N(12) = 0.25 * (1 - c(2) * c(2)) * (1 - c(0)) * (1 - c(1));
  N(13) = 0.25 * (1 - c(2) * c(2)) * (1 + c(0)) * (1 - c(1));
  N(14) = 0.25 * (1 - c(2) * c(2)) * (1 + c(0)) * (1 + c(1));
  N(15) = 0.25 * (1 - c(2) * c(2)) * (1 - c(0)) * (1 + c(1));
  N(16) = 0.25 * (1 - c(0) * c(0)) * (1 - c(1)) * (1 + c(2));
  N(17) = 0.25 * (1 - c(1) * c(1)) * (1 + c(0)) * (1 + c(2));
  N(18) = 0.25 * (1 - c(0) * c(0)) * (1 + c(1)) * (1 + c(2));
  N(19) = 0.25 * (1 - c(1) * c(1)) * (1 - c(0)) * (1 + c(2));
}
/* -------------------------------------------------------------------------- */

template <>
template <class vector_type, class matrix_type>
inline void InterpolationElement<_itp_serendip_hexahedron_20>::computeDNDS(
    const vector_type & c, matrix_type & dnds) {
  // derivatives
  // ddx
  dnds(0, 0) =
      0.25 * (c(0) + 0.5 * (c(1) + c(2) + 1)) * (c(1) - 1) * (c(2) - 1);
  dnds(0, 1) =
      0.25 * (c(0) - 0.5 * (c(1) + c(2) + 1)) * (c(1) - 1) * (c(2) - 1);
  dnds(0, 2) =
      -0.25 * (c(0) + 0.5 * (c(1) - c(2) - 1)) * (c(1) + 1) * (c(2) - 1);
  dnds(0, 3) =
      -0.25 * (c(0) - 0.5 * (c(1) - c(2) - 1)) * (c(1) + 1) * (c(2) - 1);
  dnds(0, 4) =
      -0.25 * (c(0) + 0.5 * (c(1) - c(2) + 1)) * (c(1) - 1) * (c(2) + 1);
  dnds(0, 5) =
      -0.25 * (c(0) - 0.5 * (c(1) - c(2) + 1)) * (c(1) - 1) * (c(2) + 1);
  dnds(0, 6) =
      0.25 * (c(0) + 0.5 * (c(1) + c(2) - 1)) * (c(1) + 1) * (c(2) + 1);
  dnds(0, 7) =
      0.25 * (c(0) - 0.5 * (c(1) + c(2) - 1)) * (c(1) + 1) * (c(2) + 1);
  dnds(0, 8) = -0.5 * c(0) * (c(1) - 1) * (c(2) - 1);
  dnds(0, 9) = 0.25 * (c(1) * c(1) - 1) * (c(2) - 1);
  dnds(0, 10) = 0.5 * c(0) * (c(1) + 1) * (c(2) - 1);
  dnds(0, 11) = -0.25 * (c(1) * c(1) - 1) * (c(2) - 1);
  dnds(0, 12) = -0.25 * (c(2) * c(2) - 1) * (c(1) - 1);
  dnds(0, 13) = 0.25 * (c(1) - 1) * (c(2) * c(2) - 1);
  dnds(0, 14) = -0.25 * (c(1) + 1) * (c(2) * c(2) - 1);
  dnds(0, 15) = 0.25 * (c(1) + 1) * (c(2) * c(2) - 1);
  dnds(0, 16) = 0.5 * c(0) * (c(1) - 1) * (c(2) + 1);
  dnds(0, 17) = -0.25 * (c(2) + 1) * (c(1) * c(1) - 1);
  dnds(0, 18) = -0.5 * c(0) * (c(1) + 1) * (c(2) + 1);
  dnds(0, 19) = 0.25 * (c(2) + 1) * (c(1) * c(1) - 1);

  // ddy
  dnds(1, 0) =
      0.25 * (c(1) + 0.5 * (c(0) + c(2) + 1)) * (c(0) - 1) * (c(2) - 1);
  dnds(1, 1) =
      -0.25 * (c(1) - 0.5 * (c(0) - c(2) - 1)) * (c(0) + 1) * (c(2) - 1);
  dnds(1, 2) =
      -0.25 * (c(1) + 0.5 * (c(0) - c(2) - 1)) * (c(0) + 1) * (c(2) - 1);
  dnds(1, 3) =
      0.25 * (c(1) - 0.5 * (c(0) + c(2) + 1)) * (c(0) - 1) * (c(2) - 1);
  dnds(1, 4) =
      -0.25 * (c(1) + 0.5 * (c(0) - c(2) + 1)) * (c(0) - 1) * (c(2) + 1);
  dnds(1, 5) =
      0.25 * (c(1) - 0.5 * (c(0) + c(2) - 1)) * (c(0) + 1) * (c(2) + 1);
  dnds(1, 6) =
      0.25 * (c(1) + 0.5 * (c(0) + c(2) - 1)) * (c(0) + 1) * (c(2) + 1);
  dnds(1, 7) =
      -0.25 * (c(1) - 0.5 * (c(0) - c(2) + 1)) * (c(0) - 1) * (c(2) + 1);
  dnds(1, 8) = -0.25 * (c(0) * c(0) - 1) * (c(2) - 1);
  dnds(1, 9) = 0.5 * c(1) * (c(0) + 1) * (c(2) - 1);
  dnds(1, 10) = 0.25 * (c(0) * c(0) - 1) * (c(2) - 1);
  dnds(1, 11) = -0.5 * c(1) * (c(0) - 1) * (c(2) - 1);
  dnds(1, 12) = -0.25 * (c(2) * c(2) - 1) * (c(0) - 1);
  dnds(1, 13) = 0.25 * (c(0) + 1) * (c(2) * c(2) - 1);
  dnds(1, 14) = -0.25 * (c(0) + 1) * (c(2) * c(2) - 1);
  dnds(1, 15) = 0.25 * (c(0) - 1) * (c(2) * c(2) - 1);
  dnds(1, 16) = 0.25 * (c(2) + 1) * (c(0) * c(0) - 1);
  dnds(1, 17) = -0.5 * c(1) * (c(0) + 1) * (c(2) + 1);
  dnds(1, 18) = -0.25 * (c(2) + 1) * (c(0) * c(0) - 1);
  dnds(1, 19) = 0.5 * c(1) * (c(0) - 1) * (c(2) + 1);

  // ddz
  dnds(2, 0) =
      0.25 * (c(2) + 0.5 * (c(0) + c(1) + 1)) * (c(0) - 1) * (c(1) - 1);
  dnds(2, 1) =
      -0.25 * (c(2) - 0.5 * (c(0) - c(1) - 1)) * (c(0) + 1) * (c(1) - 1);
  dnds(2, 2) =
      0.25 * (c(2) - 0.5 * (c(0) + c(1) - 1)) * (c(0) + 1) * (c(1) + 1);
  dnds(2, 3) =
      -0.25 * (c(2) + 0.5 * (c(0) - c(1) + 1)) * (c(0) - 1) * (c(1) + 1);
  dnds(2, 4) =
      0.25 * (c(2) - 0.5 * (c(0) + c(1) + 1)) * (c(0) - 1) * (c(1) - 1);
  dnds(2, 5) =
      -0.25 * (c(2) + 0.5 * (c(0) - c(1) - 1)) * (c(0) + 1) * (c(1) - 1);
  dnds(2, 6) =
      0.25 * (c(2) + 0.5 * (c(0) + c(1) - 1)) * (c(0) + 1) * (c(1) + 1);
  dnds(2, 7) =
      -0.25 * (c(2) - 0.5 * (c(0) - c(1) + 1)) * (c(0) - 1) * (c(1) + 1);
  dnds(2, 8) = -0.25 * (c(0) * c(0) - 1) * (c(1) - 1);
  dnds(2, 9) = 0.25 * (c(1) * c(1) - 1) * (c(0) + 1);
  dnds(2, 10) = 0.25 * (c(0) * c(0) - 1) * (c(1) + 1);
  dnds(2, 11) = -0.25 * (c(1) * c(1) - 1) * (c(0) - 1);
  dnds(2, 12) = -0.5 * c(2) * (c(1) - 1) * (c(0) - 1);
  dnds(2, 13) = 0.5 * c(2) * (c(0) + 1) * (c(1) - 1);
  dnds(2, 14) = -0.5 * c(2) * (c(0) + 1) * (c(1) + 1);
  dnds(2, 15) = 0.5 * c(2) * (c(0) - 1) * (c(1) + 1);
  dnds(2, 16) = 0.25 * (c(1) - 1) * (c(0) * c(0) - 1);
  dnds(2, 17) = -0.25 * (c(0) + 1) * (c(1) * c(1) - 1);
  dnds(2, 18) = -0.25 * (c(1) + 1) * (c(0) * c(0) - 1);
  dnds(2, 19) = 0.25 * (c(0) - 1) * (c(1) * c(1) - 1);
}

/* -------------------------------------------------------------------------- */

template <>
inline Real
GeometricalElement<_gt_hexahedron_20>::getInradius(const Matrix<Real> & coord) {
  Vector<Real> u0 = coord(0);
  Vector<Real> u1 = coord(1);
  Vector<Real> u2 = coord(2);
  Vector<Real> u3 = coord(3);
  Vector<Real> u4 = coord(4);
  Vector<Real> u5 = coord(5);
  Vector<Real> u6 = coord(6);
  Vector<Real> u7 = coord(7);
  Vector<Real> u8 = coord(8);
  Vector<Real> u9 = coord(9);
  Vector<Real> u10 = coord(10);
  Vector<Real> u11 = coord(11);
  Vector<Real> u12 = coord(12);
  Vector<Real> u13 = coord(13);
  Vector<Real> u14 = coord(14);
  Vector<Real> u15 = coord(15);
  Vector<Real> u16 = coord(16);
  Vector<Real> u17 = coord(17);
  Vector<Real> u18 = coord(18);
  Vector<Real> u19 = coord(19);

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

  Real x = std::min(a, std::min(b, std::min(c, d)));
  Real y = std::min(e, std::min(f, std::min(g, h)));
  Real z = std::min(i, std::min(j, std::min(k, l)));
  Real p = std::min(x, std::min(y, z));

  return p;
}
