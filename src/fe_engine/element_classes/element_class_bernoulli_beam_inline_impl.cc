/**
 * @file   element_class_bernoulli_beam_inline_impl.cc
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 *
 * @date creation: Wed Jan 16 2013
 * @date last modification: Mon Jul 07 2014
 *
 * @brief  Specialization of the element_class class for the type _bernoulli_beam_2
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
 * @section DESCRIPTION
 *
 * @verbatim
   --x-----q1----|----q2-----x---> x
    -a          0            a
 @endverbatim
 *
 * @subsection coords Nodes coordinates
 *
 * @f[
 * \begin{array}{ll}
 *  x_{1}  = -a   &  x_{2} = a
 * \end{array}
 * @f]
 *
 * @subsection shapes Shape functions
 * @f[
 *   \begin{array}{ll}
 *     N_1(x) &= \frac{1-x}{2a}\\
 *     N_2(x) &= \frac{1+x}{2a}
 *   \end{array}
 *
 *   \begin{array}{ll}
 *     M_1(x) &= 1/4(x^{3}/a^{3}-3x/a+2)\\
 *     M_2(x) &= -1/4(x^{3}/a^{3}-3x/a-2)
 *   \end{array}
 *
 *   \begin{array}{ll}
 *     L_1(x) &= a/4(x^{3}/a^{3}-x^{2}/a^{2}-x/a+1)\\
 *     L_2(x) &= a/4(x^{3}/a^{3}+x^{2}/a^{2}-x/a-1)
 *   \end{array}
 *
 *   \begin{array}{ll}
 *     M'_1(x) &= 3/4a(x^{2}/a^{2}-1)\\
 *     M'_2(x) &= -3/4a(x^{2}/a^{2}-1)
 *   \end{array}
 *
 *   \begin{array}{ll}
 *     L'_1(x) &= 1/4(3x^{2}/a^{2}-2x/a-1)\\
 *     L'_2(x) &= 1/4(3x^{2}/a^{2}+2x/a-1)
 *   \end{array}
 *@f]
 *
 * @subsection dnds Shape derivatives
 *
 *@f[
 * \begin{array}{ll}
 *   N'_1(x) &= -1/2a\\
 *   N'_2(x) &= 1/2a
 * \end{array}]
 *
 * \begin{array}{ll}
 *   -M''_1(x) &= -3x/(2a^{3})\\
 *   -M''_2(x) &= 3x/(2a^{3})\\
 * \end{array}
 *
 * \begin{array}{ll}
 *   -L''_1(x) &= -1/2a(3x/a-1)\\
 *   -L''_2(x) &= -1/2a(3x/a+1)
 * \end{array}
 *@f]
 *
 * @subsection quad_points Position of quadrature points
 *
 * @f[
 * \begin{array}{ll}
 * x_{q1}  = -a/\sqrt{3} & x_{q2} = a/\sqrt{3}
 * \end{array}
 * @f]
 */

/* -------------------------------------------------------------------------- */
AKANTU_DEFINE_STRUCTURAL_ELEMENT_CLASS_PROPERTY(_bernoulli_beam_2,
						_gt_segment_2,
						_itp_bernoulli_beam,
						_segment_2,
						_ek_structural,
						2,
						_git_segment, 4);

AKANTU_DEFINE_STRUCTURAL_ELEMENT_CLASS_PROPERTY(_bernoulli_beam_3,
						_gt_segment_2,
						_itp_bernoulli_beam,
						_segment_2,
						_ek_structural,
						3,
						_git_segment, 4);

/* -------------------------------------------------------------------------- */
template <>
inline void
InterpolationElement<_itp_bernoulli_beam>::computeShapes(const Vector<Real> & natural_coords,
							 Vector<Real> & N,
							 const Matrix<Real> & real_coord,
							 UInt id) {
  /// Compute the dimension of the beam
  Vector<Real> x1 = real_coord(0);
  Vector<Real> x2 = real_coord(1);

  Real a = .5 * x1.distance(x2);
  /// natural coordinate
  Real c = natural_coords(0);

  switch (id) {
  case 0: { // N
    N(0) = 0.5*(1 - c);
    N(1) = 0.5*(1 + c);
    break;
  }
  case 1: { // M
    N(0) =  0.25 * (c*c*c - 3*c + 2);
    N(1) = -0.25 * (c*c*c - 3*c - 2);
    break;
  }
  case 2: { // L
    N(0) = 0.25*a * (c*c*c - c*c - c + 1);
    N(1) = 0.25*a * (c*c*c + c*c - c - 1);
    break;
  }
  case 3: { // M'
    N(0) =  0.75/a * (c*c - 1);
    N(1) = -0.75/a * (c*c - 1);
    break;
  }
  case 4: { // L'
    N(0) = 0.25 * (3*c*c - 2*c - 1);
    N(1) = 0.25 * (3*c*c + 2*c - 1);
    break;
  }
  }
 }

/* -------------------------------------------------------------------------- */
template <>
inline void
InterpolationElement<_itp_bernoulli_beam>::computeDNDS(const Vector<Real> & natural_coords,
						       Matrix<Real> & dnds,
						       const Matrix<Real> & real_nodes_coord,
						       UInt id) {
  /// Compute the dimension of the beam
  Vector<Real> x1 = real_nodes_coord(0);
  Vector<Real> x2 = real_nodes_coord(1);

  Real a = .5 * x1.distance(x2);
  /// natural coordinate
  Real c = natural_coords(0)*a;


  switch (id) {
  case 0: { // N'
    dnds(0, 0) = -0.5/a;
    dnds(0, 1) =  0.5/a;
    break;
  }
  case 1: { // M''
    dnds(0, 0) = -3.*c/(2.*pow(a,3));
    dnds(0, 1) =  3.*c/(2.*pow(a,3));
    break;
  }
  case 2: { // L''
    dnds(0, 0) = -0.5/a * (3*c/a - 1);
    dnds(0, 1) =- 0.5/a * (3*c/a + 1);
    break;
  }
  }
}
