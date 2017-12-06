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
//#include "aka_element_classes_info.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_ELEMENT_CLASS_BERNOULLI_BEAM_INLINE_IMPL_CC__
#define __AKANTU_ELEMENT_CLASS_BERNOULLI_BEAM_INLINE_IMPL_CC__

namespace akantu {
/* -------------------------------------------------------------------------- */
AKANTU_DEFINE_STRUCTURAL_INTERPOLATION_TYPE_PROPERTY(_itp_bernoulli_beam_2,
                                                     _itp_lagrange_segment_2,
                                                     3, 2);

AKANTU_DEFINE_STRUCTURAL_INTERPOLATION_TYPE_PROPERTY(_itp_bernoulli_beam_3,
                                                     _itp_lagrange_segment_2,
                                                     6, 4);

AKANTU_DEFINE_STRUCTURAL_ELEMENT_CLASS_PROPERTY(_bernoulli_beam_2,
                                                _gt_segment_2,
                                                _itp_bernoulli_beam_2,
                                                _segment_2, _ek_structural, 2,
                                                _git_segment, 3);

AKANTU_DEFINE_STRUCTURAL_ELEMENT_CLASS_PROPERTY(_bernoulli_beam_3,
                                                _gt_segment_2,
                                                _itp_bernoulli_beam_3,
                                                _segment_2, _ek_structural, 3,
                                                _git_segment, 3);

/* -------------------------------------------------------------------------- */
template <>
inline void
InterpolationElement<_itp_bernoulli_beam_2, _itk_structural>::computeShapes(
    const Vector<Real> & natural_coords, Matrix<Real> & N) {
  Vector<Real> L(2);
  InterpolationElement<_itp_lagrange_segment_2, _itk_lagrangian>::computeShapes(
      natural_coords, L);
  Matrix<Real> H(2, 4);
  InterpolationElement<_itp_hermite_2, _itk_structural>::computeShapes(
    natural_coords, H);

  // clang-format off
  //    u1   v1      t1      u2   v2      t2
  N = {{L(0), 0      , 0      , L(1), 0      , 0      },  // u
       {0   , H(0, 0), H(0, 1), 0   , H(0, 2), H(0, 3)},  // v
       {0   , H(1, 0), H(1, 1), 0   , H(1, 2), H(1, 3)}}; // theta
  // clang-format on
}

template <>
inline void
InterpolationElement<_itp_bernoulli_beam_3, _itk_structural>::computeShapes(
    const Vector<Real> & natural_coords, Matrix<Real> & N) {
  Vector<Real> L(2);
  InterpolationElement<_itp_lagrange_segment_2, _itk_lagrangian>::computeShapes(
      natural_coords, L);
  Matrix<Real> H(2, 4);
  InterpolationElement<_itp_hermite_2, _itk_structural>::computeShapes(
      natural_coords, H);

  // clang-format off
  //   u1    v1      w1      x1   y1      z1      u2   v2      w2      x2   y2      z2
  N = {{L(0), 0      , 0      , 0   , 0      , 0      , L(1), 0      , 0      , 0   , 0      , 0      },  // u
       {0   , H(0, 0), 0      , 0   , H(0, 1), 0      , 0   , H(0, 2), 0      , 0   , H(0, 3), 0      },  // v
       {0   , 0      , H(0, 0), 0   , 0      , H(0, 1), 0   , 0      , H(0, 2), 0   , 0      , H(0, 3)},  // w
       {0   , 0      , 0      , L(0), 0      , 0      , 0   , 0      , 0      , L(1), 0      , 0      },  // thetax
       {0   , H(1, 0), 0      , 0   , H(1, 1), 0      , 0   , H(1, 2), 0      , 0   , H(1, 3), 0      },  // thetay
       {0   , 0      , H(1, 0), 0   , 0      , H(1, 1), 0   , 0      , H(1, 2), 0   , 0      , H(1, 3)}}; // thetaz
  // clang-format on
}


/* -------------------------------------------------------------------------- */
template <>
inline void
InterpolationElement<_itp_bernoulli_beam_2, _itk_structural>::computeDNDS(
    const Vector<Real> & natural_coords, Matrix<Real> & B) {
  Matrix<Real> L(1, 2);
  InterpolationElement<_itp_lagrange_segment_2, _itk_lagrangian>::computeDNDS(
      natural_coords, L);
  Matrix<Real> H(1, 4);
  InterpolationElement<_itp_hermite_2, _itk_structural>::computeDNDS(
      natural_coords, H);

  // clang-format off
  //    u1      v1      t1      u2      v2      t2
  B = {{L(0, 0), 0      , 0      , L(0, 1), 0      , 0      },   // epsilon
       {0      , H(0, 0), H(0, 1), 0      , H(0, 2), H(0, 3)}};  // chi (curv.)
  // clang-format on
}

template <>
inline void
InterpolationElement<_itp_bernoulli_beam_3, _itk_structural>::computeDNDS(
    const Vector<Real> & natural_coords, Matrix<Real> & B) {
  Matrix<Real> L(1, 2);
  InterpolationElement<_itp_lagrange_segment_2, _itk_lagrangian>::computeDNDS(
      natural_coords, L);
  Matrix<Real> H(1, 4);
  InterpolationElement<_itp_hermite_2, _itk_structural>::computeDNDS(
      natural_coords, H);
  // clang-format off
  //    u1      v1      w1      x1      y1      z1      u2      v2      w2      x2      y2      z2
  B = {{L(0, 0), 0      , 0      , 0      , 0      , 0      , L(0, 1), 0      , 0      , 0      , 0      , 0      },  // eps
       {0      , H(0, 0), 0      , 0      , H(0, 1), 0      , 0      , H(0, 2), 0      , 0      , H(0, 3), 0      },  // chix
       {0      , 0      , H(0, 0), 0      , 0      , H(0, 1), 0      , 0      , H(0, 2), 0      , 0      , H(0, 3)},  // chiy
       {0      , 0      , 0      , L(0, 0), 0      , 0      , 0      , 0      , 0      , L(0, 1), 0      , 0      }}; // chiz
  // clang-format on
}

template <>
inline void ElementClass<_bernoulli_beam_2>::computeRotation(
    const Matrix<Real> & node_coords, Matrix<Real> & rotation) {
  auto X1 = node_coords(0);
  auto X2 = node_coords(1);
  auto vec = Vector<Real>(X2) - Vector<Real>(X1);
  auto L = vec.norm();
  auto c = vec(0) / L;
  auto s = vec(1) / L;

  rotation = {{c, -s}, {s, c}};
}

} // namespace akantu
#endif /* __AKANTU_ELEMENT_CLASS_BERNOULLI_BEAM_INLINE_IMPL_CC__ */
