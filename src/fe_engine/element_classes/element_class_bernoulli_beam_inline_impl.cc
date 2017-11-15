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
    -1          0            1
 @endverbatim
 *
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
    void computeShapes(const Vector<Real> & natural_coords, Matrix<Real> & N) {
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
    void computeDNDS(const Vector<Real> & natural_coords, Matrix<Real> & B) {

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
    const Vector<Real> & natural_coords, Matrix<Real> & N) {
  Vector<Real> L(2);
  InterpolationElement<_itp_segment_2, _itk_lagrange>::computeShapes(
      natural_coords, L);
  Matrix<Real> H(2, 4);
  InterpolationElement<_itp_hermite_2, _itk_structural>::computeShapes(
      natural_coords, H);

  // clang-format off
  //    u1   v1      t1      u2   v2      t2
  N = {{L(0) 0       0       L(1) 0       0      },  // u
       {0    H(0, 0) H(0, 1) 0    H(0, 2) H(0, 3)},  // v
       {0    H(1, 0) H(1, 1) 0    H(1, 2) H(1, 3)}}; // theta
  // clang-format on
}

template <>
inline void
InterpolationElement<_itp_bernoulli_beam_3, _itk_structural>::computeShapes(
    const Vector<Real> & natural_coords, Matrix<Real> & N) {
  Vector<Real> L(2);
  InterpolationElement<_itp_segment_2, _itk_lagrange>::computeShapes(
      natural_coords, L);
  Matrix<Real> H(2, 4);
  InterpolationElement<_itp_hermite_2, _itk_structural>::computeShapes(
      natural_coords, H);

  // clang-format off
  //   u1    v1      w1      x1   y1      z1      u2   v2      w2      x2   y2      z2
  N = {{L(0) 0       0       0    0       0       L(1) 0       0       0    0       0      },  // u
       {0    H(0, 0) 0       0    H(0, 1) 0       0    H(0, 2) 0       0    H(0, 3) 0      },  // v
       {0    0       H(0, 0) 0    0       H(0, 1) 0    0       H(0, 2) 0    0       H(0, 3)},  // w
       {0    0       0       L(0) 0       0       0    0       0       L(1) 0       0      },  // thetax
       {0    H(1, 0) 0       0    H(1, 1) 0       0    H(1, 2) 0       0    H(1, 3) 0      },  // thetay
       {0    0       H(1, 0) 0    0       H(1, 1) 0    0       H(1, 2) 0    0       H(1, 3)}}; // thetaz
  // clang-format on
}


/* -------------------------------------------------------------------------- */
template <>
inline void
InterpolationElement<_itp_bernoulli_beam_2, _itk_structural>::computeDNDS(
    const Vector<Real> & natural_coords, Matrix<Real> & B) {
  Matrix<Real> L(1, 2);
  InterpolationElement<_itp_segment_2, _itk_lagrange>::computeDNDS(
      natural_coords, L);
  Matrix<Real> H(1, 4);
  InterpolationElement<_itp_hermite_2, _itk_structural>::computeDNDS(
      natural_coords, H);

  // clang-format off
  //    u1      v1      t1      u2      v2      t2
  B = {{L(0, 0) 0       0       L(0, 1) 0       0      },   // epsilon
       {0       H(0, 0) H(0, 1) 0       H(0, 2) H(0, 3)}};  // chi (curv.)
  // clang-format on
}

template <>
inline void
InterpolationElement<_itp_bernoulli_beam_3, _itk_structural>::computeDNDS(
    const Vector<Real> & natural_coords, Matrix<Real> & B) {
  Matrix<Real> L(1, 2);
  InterpolationElement<_itp_segment_2, _itk_lagrange>::computeDNDS(
      natural_coords, L);
  Matrix<Real> H(1, 4);
  InterpolationElement<_itp_hermite_2, _itk_structural>::computeDNDS(
      natural_coords, H);
  // clang-format off
  //    u1      v1      w1      x1      y1      z1      u2      v2      w2      x2      y2      z2
  B = {{L(0, 0) 0       0       0       0       0       L(0, 1) 0       0       0       0       0      },  // eps
       {0       H(0, 0) 0       0       H(0, 1) 0       0       H(0, 2) 0       0       H(0, 3) 0      },  // chix
       {0       0       H(0, 0) 0       0       H(0, 1) 0       0       H(0, 2) 0       0       H(0, 3)},  // chiy
       {0       0       0       L(0, 0) 0       0       0       0       0       L(0, 1) 0       0      }}; // chiz
  // clang-format on
}

} // namespace akantu
#endif /* __AKANTU_ELEMENT_CLASS_BERNOULLI_BEAM_INLINE_IMPL_CC__ */
