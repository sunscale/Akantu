/**
 * @file   element_class_pentahedron_6_inline_impl.cc
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
 * @section LICENSE
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
 * @section DESCRIPTION
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
 */

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
  dnds(0, 3) = 0.5 * c(1);
  dnds(0, 4) = 0.5 * c(2);
  dnds(0, 5) = 0.5 * (1 - c(1) - c(2));

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
template <>
inline Real
GeometricalElement<_gt_pentahedron_6>::getInradius(const Matrix<Real> & coord) {
  Vector<Real> u0 = coord(0);
  Vector<Real> u1 = coord(1);
  Vector<Real> u2 = coord(2);
  Vector<Real> u3 = coord(3);

  Real a = u0.distance(u1);
  Real b = u1.distance(u2);
  Real c = u2.distance(u3);
  Real d = u3.distance(u0);
  Real s = (a + b + c) / 2;
  Real A = std::sqrt(s * (s - a) * (s - b) * (s - c));
  Real ra = 2 * s / A;
  Real p = std::min(ra, d);

  return p;
}
