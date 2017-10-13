/**
 * @file   element_class_bernoulli_beam_inline_impl.cc
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 *
 * @date creation: Fri Jul 15 2011
 * @date last modification: Sun Oct 19 2014
 *
 * @brief  Specialization of the element_class class for the type
 _bernoulli_beam_2
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
#include "aka_static_if.hh"
#include "element_class_structural.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_ELEMENT_CLASS_BERNOULLI_BEAM_INLINE_IMPL_CC__
#define __AKANTU_ELEMENT_CLASS_BERNOULLI_BEAM_INLINE_IMPL_CC__

namespace akantu {
/* -------------------------------------------------------------------------- */
AKANTU_DEFINE_STRUCTURAL_INTERPOLATION_TYPE_PROPERTY(_itp_bernoulli_beam_2,
                                                     _itp_lagrange_segment_2, 2,
                                                     3, 2);
AKANTU_DEFINE_STRUCTURAL_INTERPOLATION_TYPE_PROPERTY(_itp_bernoulli_beam_3,
                                                     _itp_lagrange_segment_2, 3,
                                                     6, 4);

AKANTU_DEFINE_STRUCTURAL_ELEMENT_CLASS_PROPERTY(_bernoulli_beam_2,
                                                _gt_segment_2,
                                                _itp_bernoulli_beam_2,
                                                _segment_2, _ek_structural, 2,
                                                _git_segment, 5);

AKANTU_DEFINE_STRUCTURAL_ELEMENT_CLASS_PROPERTY(_bernoulli_beam_3,
                                                _gt_segment_2,
                                                _itp_bernoulli_beam_3,
                                                _segment_2, _ek_structural, 3,
                                                _git_segment, 5);

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
namespace {
  namespace details {
    template <InterpolationType type>
    void computeShapes(const Vector<Real> & natural_coords, Matrix<Real> & N,
                       const Matrix<Real> & real_coord) {
      /// Compute the dimension of the beam
      Vector<Real> x1 = real_coord(0);
      Vector<Real> x2 = real_coord(1);

      Real a = x1.distance(x2) / 2.;
      /// natural coordinate
      Real c = natural_coords(0);

      auto N0 = (1 - c) / 2.;
      auto N1 = (1 + c) / 2.;

      auto M0 = (c * c * c - 3. * c + 2.) / 4.;
      auto M1 = -(c * c * c - 3. * c - 2.) / 4.;

      auto L0 = a * (c * c * c - c * c - c + 1.) / 4.;
      auto L1 = a * (c * c * c + c * c - c - 1.) / 4.;

      auto Mp0 = 3. / a * (c * c - 1.) / 4.;
      auto Mp1 = -3. / a * (c * c - 1.) / 4.;

      auto Lp0 = (3. * c * c - 2. * c - 1.) / 4.;
      auto Lp1 = (3. * c * c + 2. * c - 1.) / 4.;

      static_if(type == _itp_bernoulli_beam_2)
          .then([&](auto && N) {
            // clang-format off
            //     0    1    2   3    4    5
            N = {{N0,  0.,  0., N1,  0.,  0.},
                 {0.,  M0,  L0, 0.,  M1,  L1},
                 {0., Mp0, Lp0, 0., Mp1, Lp1}};
            // clang-format on
          })
          .else_if(type == _itp_bernoulli_beam_3)
          .then([&](auto && N) {
            // clang-format off
            //     0    1    2    3     4    5    6    7    8   9    10   11
            N = {{N0,   0.,  0.,  0.,   0.,  0., N1,   0.,  0., 0.,   0.,  0.},
                 { 0., M0,   0.,  0.,   0.,  L0,  0.,  M1,  0., 0.,   0.,  L1},
                 { 0.,  0.,  M0,  0.,  -L0,  0.,  0.,  0.,  M1, 0.,  -L1,  0.},
                 { 0.,  0.,  0., N0,    0.,  0.,  0.,  0.,  0., N1,   0.,  0.},
                 { 0.,  0., Mp0,  0., -Lp0,  0.,  0.,  0., Mp1, 0., -Lp1,  0.},
                 { 0., Mp0,  0.,  0.,   0., Lp0,  0., Mp1,  0., 0.,   0., Lp1}};
            // clang-format on
          })
          .else_([](auto && /*unused*/) {
            AKANTU_EXCEPTION("Should not be in this part of the code");
          })(std::forward<decltype(N)>(N));
    }

    /* ---------------------------------------------------------------------- */

    template <InterpolationType type>
    void computeDNDS(const Vector<Real> & natural_coords, Matrix<Real> & B,
                     const Matrix<Real> & real_nodes_coord) {
      /// Compute the dimension of the beam
      Vector<Real> x1 = real_nodes_coord(0);
      Vector<Real> x2 = real_nodes_coord(1);

      Real a = .5 * x1.distance(x2);
      /// natural coordinate
      Real c = natural_coords(0) * a;

      auto Np0 = -1. / (2. * a);
      auto Np1 = 1. / (2. * a);

      auto Mpp0 = -3. * c / (2. * pow(a, 3.));
      auto Mpp1 = 3. * c / (2. * pow(a, 3.));

      auto Lpp0 = -1. / (2. * a) * (3. * c / a - 1.);
      auto Lpp1 = -1. / (2. * a) * (3. * c / a + 1.);

      static_if(type == _itp_bernoulli_beam_2)
          .then([&](auto && B) {
            // clang-format off
            //      0    1     2    3     4     5
            B = {{Np0,   0.,   0., Np1,   0.,   0.},
                 { 0., Mpp0, Lpp0,  0., Mpp1, Lpp1}};
            // clang-format on
          })
          .else_if(type == _itp_bernoulli_beam_3)
          .then([&](auto && B) {
            // clang-format off
            //      0    1     2    3      4     5     6    7     8    9     10    11
            B = {{Np0,   0.,   0.,  0.,    0.,   0., Np1,   0.,   0.,  0.,    0.,   0.},
                 { 0., Mpp0,   0.,  0.,    0., Lpp0,  0., Mpp1,   0.,  0.,    0., Lpp1},
                 { 0.,   0., Mpp0,  0., -Lpp0,   0.,  0.,   0., Mpp1,  0., -Lpp1,   0.},
                 { 0.,   0.,   0., Np0,    0.,   0.,  0.,   0.,   0., Np1,    0.,   0.}};
            // clang-format on
          })
          .else_([](auto && /*unused*/) {
            AKANTU_EXCEPTION("Should not be in this part of the code");
          })(std::forward<decltype(B)>(B));
    }
  } // namespace details
} // namespace

/* -------------------------------------------------------------------------- */
template <>
inline void
InterpolationElement<_itp_bernoulli_beam_2, _itk_structural>::computeShapes(
    const Vector<Real> & natural_coords, Matrix<Real> & N,
    const Matrix<Real> & real_coord) {
  details::computeShapes<_itp_bernoulli_beam_2>(natural_coords, N, real_coord);
}

template <>
inline void
InterpolationElement<_itp_bernoulli_beam_3, _itk_structural>::computeShapes(
    const Vector<Real> & natural_coords, Matrix<Real> & N,
    const Matrix<Real> & real_coord) {
  details::computeShapes<_itp_bernoulli_beam_3>(natural_coords, N, real_coord);
}

/* -------------------------------------------------------------------------- */
template <>
inline void
InterpolationElement<_itp_bernoulli_beam_2, _itk_structural>::computeDNDS(
    const Vector<Real> & natural_coords, Matrix<Real> & B,
    const Matrix<Real> & real_nodes_coord) {
  details::computeDNDS<_itp_bernoulli_beam_2>(natural_coords, B,
                                              real_nodes_coord);
}

template <>
inline void
InterpolationElement<_itp_bernoulli_beam_3, _itk_structural>::computeDNDS(
    const Vector<Real> & natural_coords, Matrix<Real> & B,
    const Matrix<Real> & real_nodes_coord) {
  details::computeDNDS<_itp_bernoulli_beam_3>(natural_coords, B,
                       real_nodes_coord);
}

} // namespace akantu
#endif /* __AKANTU_ELEMENT_CLASS_BERNOULLI_BEAM_INLINE_IMPL_CC__ */
