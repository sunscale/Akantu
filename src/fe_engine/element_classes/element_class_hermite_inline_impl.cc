/**
 * @file   element_class_hermite_inline_impl.cc
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Fri Jul 15 2011
 * @date last modification: Sun Oct 19 2014
 *
 * @brief  Specialization of the element_class class for the type
 _hermite
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
    -1          0            1
 @endverbatim
 *
 * @subsection shapes Shape functions
 * @f[
 *   \begin{array}{ll}
 *     M_1(\xi) &= 1/4(\xi^{3}/-3\xi+2)\\
 *     M_2(\xi) &= -1/4(\xi^{3}-3\xi-2)
 *   \end{array}
 *
 *   \begin{array}{ll}
 *     L_1(\xi) &= 1/4(\xi^{3}-\xi^{2}-\xi+1)\\
 *     L_2(\xi) &= 1/4(\xi^{3}+\xi^{2}-\xi-1)
 *   \end{array}
 *
 *   \begin{array}{ll}
 *     M'_1(\xi) &= 3/4(\xi^{2}-1)\\
 *     M'_2(\xi) &= -3/4(\xi^{2}-1)
 *   \end{array}
 *
 *   \begin{array}{ll}
 *     L'_1(\xi) &= 1/4(3\xi^{2}-2\xi-1)\\
 *     L'_2(\xi) &= 1/4(3\xi^{2}+2\xi-1)
 *   \end{array}
 *@f]
 *
 * @subsection dnds Shape derivatives
 *
 *@f[
 * \begin{array}{ll}
 *   N'_1(\xi) &= -1/2\\
 *   N'_2(\xi) &= 1/2
 * \end{array}]
 *
 * \begin{array}{ll}
 *   -M''_1(\xi) &= -3\xi/2\\
 *   -M''_2(\xi) &= 3\xi/2\\
 * \end{array}
 *
 * \begin{array}{ll}
 *   -L''_1(\xi) &= -1/2a(3\xi/a-1)\\
 *   -L''_2(\xi) &= -1/2a(3\xi/a+1)
 * \end{array}
 *@f]
 *
 */
/* -------------------------------------------------------------------------- */
#include "aka_static_if.hh"
#include "element_class_structural.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_ELEMENT_CLASS_HERMITE_INLINE_IMPL_CC__
#define __AKANTU_ELEMENT_CLASS_HERMITE_INLINE_IMPL_CC__

namespace akantu {
/* -------------------------------------------------------------------------- */

AKANTU_DEFINE_STRUCTURAL_INTERPOLATION_TYPE_PROPERTY(_itp_hermite_2,
                                                     _itp_lagrange_segment_2, 2,
                                                     1);

/* -------------------------------------------------------------------------- */

namespace {
  namespace details {
    template <InterpolationType type>
    void computeShapes(const Vector<Real> & natural_coords, Matrix<Real> & N) {
      /// natural coordinate
      Real xi = natural_coords(0);

      // Cubic Hermite splines interpolating displacement
      auto M1 = 1. / 4. * Math::pow<2>(xi - 1) * (xi + 2);
      auto M2 = 1. / 4. * Math::pow<2>(xi + 1) * (2 - xi);
      auto L1 = 1. / 4. * Math::pow<2>(xi - 1) * (xi + 1);
      auto L2 = 1. / 4. * Math::pow<2>(xi + 1) * (xi - 1);

#if 1 // Version where we also interpolate the rotations
      // Derivatives of previous functions interpolating rotations
      auto M1_ = 3. / 4. * (xi * xi - 1);
      auto L1_ = 1. / 4. * (3 * xi * xi - 2 * xi - 1);
      auto M2_ = 3. / 4. * (1 - xi * xi);
      auto L2_ = 1. / 4. * (3 * xi * xi + 2 * xi - 1);

      // clang-format off
      //    v1   t1   v2   t2
      N = {{M1 , L1,  M2,  L2}, // displacement interpolation
	   {M1_, L1_, M2_, L2_}}; // rotation interpolation
      // clang-format on

#else // Version where we only interpolate displacements
      // clang-format off
      //    v1  t1  v2  t2
      N = {{M1, L1, M2, L2}};
      // clang-format on
#endif
    }

    /* ---------------------------------------------------------------------- */

    template <InterpolationType type>
    void computeDNDS(const Vector<Real> & natural_coords, Matrix<Real> & B) {
      // natural coordinate
      Real xi = natural_coords(0);
      // Derivatives for rotations
      auto M1__ = 3. / 2. * xi;
      auto L1__ = 1. / 2. * (3 * xi - 1);
      auto M2__ = 3. / 2. * (-xi);
      auto L2__ = 1. / 2. * (3 * xi + 1);

      //    v1    t1    v2    t2
      B = {{M1__, L1__, M2__, L2__}}; // computing curvature : {chi} = [B]{d}
    }
  } // namespace details
} // namespace

/* -------------------------------------------------------------------------- */
template <>
inline void
InterpolationElement<_itp_hermite_2, _itk_structural>::computeShapes(
    const Vector<Real> & natural_coords, Matrix<Real> & N) {
  details::computeShapes<_itp_hermite_2>(natural_coords, N);
}

/* -------------------------------------------------------------------------- */
template <>
inline void
InterpolationElement<_itp_hermite_2, _itk_structural>::computeDNDS(
    const Vector<Real> & natural_coords, Matrix<Real> & B) {
  details::computeDNDS<_itp_hermite_2>(natural_coords, B);
}

} // namespace akantu
#endif /* __AKANTU_ELEMENT_CLASS_HERMITE_INLINE_IMPL_CC__ */
