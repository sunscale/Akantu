/**
 * @file   element_class_pentahedron_6_inline_impl.cc
 *
 * @author Marion Estelle Chambart <mchambart@stucky.ch>
 * @author Thomas Menouillard <tmenouillard@stucky.ch>
 *
 * @date creation: Wed Jun 12 2013
 * @date last modification: Fri Jun 13 2014
 *
 * @brief  Specialization of the element_class class for the type _pentahedron_6
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
                   \zeta
                    ^
         (-1,1,1)   |     (1,1,1)
                8---|------7
               /|   |     /|
              / |   |    / |
   (-1,-1,1) 5----------6  | (1,-1,1)
             |  |   |   |  |
             |  |   |   |  |
             |  |   +---|-------> \xi
             |  |  /    |  |
   (-1,1,-1) |  4-/-----|--3 (1,1,-1)
             | / /      | /
             |/ /       |/
             1-/--------2
   (-1,-1,-1) /        (1,-1,-1)
             /
            \eta
 @endverbatim
 *
 * @subsection shapes Shape functions
 * @f[
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
 * N43 = (1 - \xi) (1 + \eta) (1 - \zeta) / 8
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
 * @f]
 *
 * @subsection quad_points Position of quadrature points
 * @f{eqnarray*}{
 * \xi_{q0}  &=& -1/\sqrt{3} \qquad  \eta_{q0} = -1/\sqrt{3} \qquad \zeta_{q0} = -1/\sqrt{3} \\
 * \xi_{q1}  &=&  1/\sqrt{3} \qquad  \eta_{q1} = -1/\sqrt{3} \qquad \zeta_{q1} = -1/\sqrt{3} \\
 * \xi_{q2}  &=&  1/\sqrt{3} \qquad  \eta_{q2} =  1/\sqrt{3} \qquad \zeta_{q2} = -1/\sqrt{3} \\
 * \xi_{q3}  &=& -1/\sqrt{3} \qquad  \eta_{q3} =  1/\sqrt{3} \qquad \zeta_{q3} = -1/\sqrt{3} \\
 * \xi_{q4}  &=& -1/\sqrt{3} \qquad  \eta_{q4} = -1/\sqrt{3} \qquad \zeta_{q4} =  1/\sqrt{3} \\
 * \xi_{q5}  &=&  1/\sqrt{3} \qquad  \eta_{q5} = -1/\sqrt{3} \qquad \zeta_{q5} =  1/\sqrt{3} \\
 * \xi_{q6}  &=&  1/\sqrt{3} \qquad  \eta_{q6} =  1/\sqrt{3} \qquad \zeta_{q6} =  1/\sqrt{3} \\
 * \xi_{q7}  &=& -1/\sqrt{3} \qquad  \eta_{q7} =  1/\sqrt{3} \qquad \zeta_{q7} =  1/\sqrt{3} \\
 * @f}
 */

/* -------------------------------------------------------------------------- */
AKANTU_DEFINE_ELEMENT_CLASS_PROPERTY(_pentahedron_6,
				     _gt_pentahedron_6,
				     _itp_lagrange_pentahedron_6,
				     _ek_regular,
				     3,
				     _git_pentahedron,
				     1);

AKANTU_DEFINE_SHAPE(_gt_pentahedron_6, _gst_prism);

/* -------------------------------------------------------------------------- */
template <>
template <class vector_type>
inline void
InterpolationElement<_itp_lagrange_pentahedron_6>::computeShapes(const vector_type & c,
								 vector_type & N) {
  /// Natural coordinates
  N(0) =  0.5*c(0)*(1-c(2));           // N1(q)
  N(1) =  0.5*c(1)*(1-c(2));           // N2(q)
  N(2) =  0.5*(1-c(0)-c(1))*(1-c(2));  // N3(q)
  N(3) =  0.5*c(0)*(c(2)+1);           // N4(q)
  N(4) =  0.5*c(1)*(c(2)+1);           // N5(q)
  N(5) =  0.5*(1-c(0)-c(1))*(c(2)+1);  // N6(q)
}
/* -------------------------------------------------------------------------- */
template <>
template <class vector_type, class matrix_type>
inline void
InterpolationElement<_itp_lagrange_pentahedron_6>::computeDNDS(const vector_type & c,
                                                              matrix_type & dnds) {
  /**
   * @f[
   * dnds = \left(
   *          \begin{array}{cccccccc}
   *            \frac{\partial N1}{\partial \xi}  & \frac{\partial N2}{\partial \xi}
   *               & \frac{\partial N3}{\partial \xi}  & \frac{\partial N4}{\partial \xi}
   *               & \frac{\partial N5}{\partial \xi}  & \frac{\partial N6}{\partial \xi}
   *               & \frac{\partial N7}{\partial \xi}  & \frac{\partial N8}{\partial \xi}\\
   *            \frac{\partial N1}{\partial \eta} & \frac{\partial N2}{\partial \eta}
   *               & \frac{\partial N3}{\partial \eta} & \frac{\partial N4}{\partial \eta}
   *               & \frac{\partial N5}{\partial \eta} & \frac{\partial N6}{\partial \eta}
   *               & \frac{\partial N7}{\partial \eta} & \frac{\partial N8}{\partial \eta}\\
   *            \frac{\partial N1}{\partial \zeta} & \frac{\partial N2}{\partial \zeta}
   *               & \frac{\partial N3}{\partial \zeta} & \frac{\partial N4}{\partial \zeta}
   *               & \frac{\partial N5}{\partial \zeta} & \frac{\partial N6}{\partial \zeta}
   *               & \frac{\partial N7}{\partial \zeta} & \frac{\partial N8}{\partial \zeta}
   *          \end{array}
   *        \right)
   * @f]
   */
  dnds(0, 0) = 0.5*(1-c(2));
  dnds(0, 1) = 0 ;
  dnds(0, 2) = -0.5*(1-c(2));
  dnds(0, 3) =  0.5*(c(2)+1);
  dnds(0, 4) =  0.;
  dnds(0, 5) =  -0.5*(1+c(2));

  dnds(1, 0) =  0. ;
  dnds(1, 1) =  0.5*(1-c(2));
  dnds(1, 2) = -0.5*(1-c(2));
  dnds(1, 3) =  0.;
  dnds(1, 4) =  0.5*(c(2)+1);
  dnds(1, 5) = -0.5*(1+c(2));
 
  dnds(2, 0) =  -0.5*c(0);
  dnds(2, 1) =  -0.5*c(1);
  dnds(2, 2) = -0.5*(1-c(0)-c(1));
  dnds(2, 3) =  0.5*c(0);
  dnds(2, 4) =  0.5*c(1);
  dnds(2, 5) = 0.5*(1-c(0)-c(1));
}

/* -------------------------------------------------------------------------- */
template<>
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
  Real s = (a+b+c)/2;
  Real A = std::sqrt(s*(s-a)*(s-b)*(s-c));
  Real ra = 2*s/A;
  Real p = std::min(ra, d);

  return p;
}
