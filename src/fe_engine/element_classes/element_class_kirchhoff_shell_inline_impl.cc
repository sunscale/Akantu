/**
 * @file   element_class_kirchhoff_shell_inline_impl.cc
 *
 * @author Damien Spielmann <damien.spielmann@epfl.ch>
 *
 * @date creation: Fri Jul 04 2014
 * @date last modification: Sun Oct 19 2014
 *
 * @brief  Element class Kirchhoff Shell
 *
 * @section LICENSE
 *
 * Copyright  (©)  2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de Lausanne)
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
 */
/* -------------------------------------------------------------------------- */
#include "element_class_structural.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_ELEMENT_CLASS_KIRCHHOFF_SHELL_INLINE_IMPL_CC__
#define __AKANTU_ELEMENT_CLASS_KIRCHHOFF_SHELL_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
AKANTU_DEFINE_STRUCTURAL_INTERPOLATION_TYPE_PROPERTY(_itp_discrete_kirchhoff_triangle_18,
                                                     _itp_lagrange_triangle_3,
                                                     6, 6);
AKANTU_DEFINE_STRUCTURAL_ELEMENT_CLASS_PROPERTY(_discrete_kirchhoff_triangle_18,
                                                _gt_triangle_3,
                                                _itp_discrete_kirchhoff_triangle_18,
                                                _triangle_3, _ek_structural, 3,
                                                _git_triangle, 2);

/* -------------------------------------------------------------------------- */
template <>
inline void InterpolationElement<_itp_discrete_kirchhoff_triangle_18>::computeShapes(
    const Vector<Real> & /*natural_coords*/, Matrix<Real> & /*N*/) {
  //   // projected_coord (x1 x2 x3) (y1 y2 y3)

#if 0
  // natural coordinate
  Real xi = natural_coords(0);
  Real eta = natural_coords(1);

  Real x21 = projected_coord(0, 0) - projected_coord(0, 1); // x1-x2
  Real x32 = projected_coord(0, 1) - projected_coord(0, 2);
  Real x13 = projected_coord(0, 2) - projected_coord(0, 0);

  Real y21 = projected_coord(1, 0) - projected_coord(1, 1); // y1-y2
  Real y32 = projected_coord(1, 1) - projected_coord(1, 2);
  Real y13 = projected_coord(1, 2) - projected_coord(1, 0);

  /* Real x21=projected_coord(0,1)-projected_coord(0,0);
  Real x32=projected_coord(0,2)-projected_coord(0,1);
  Real x13=projected_coord(0,0)-projected_coord(0,2);

  Real y21=projected_coord(1,1)-projected_coord(1,0);
  Real y32=projected_coord(1,2)-projected_coord(1,1);
  Real y13=projected_coord(1,0)-projected_coord(1,2);*/

  // natural triangle side length
  Real L4 = std::sqrt(x21 * x21 + y21 * y21);
  Real L5 = std::sqrt(x32 * x32 + y32 * y32);
  Real L6 = std::sqrt(x13 * x13 + y13 * y13);

  // sinus and cosinus
  Real C4 = x21 / L4; // 1.
  Real C5 = x32 / L5; //-1./std::sqrt(2);
  Real C6 = x13 / L6; // 0.
  Real S4 = y21 / L4; // 0.;
  Real S5 = y32 / L5; // 1./std::sqrt(2);
  Real S6 = y13 / L6; //-1.;

  Real N1 = 1. - xi - eta;
  Real N2 = xi;
  Real N3 = eta;

  Real P4 = 4. * xi * (1. - xi - eta);
  Real P5 = 4. * xi * eta;
  Real P6 = 4. * eta * (1. - xi - eta);

  Vector<Real> N0{N1, N2, N3};
  Vector<Real> Nw2{-(1. / 8.) * P4 * L4 * C4 + (1. / 8.) * P6 * L6 * C6,
                   -(1. / 8.) * P5 * L5 * C5 + (1. / 8.) * P4 * L4 * C4,
                   -(1. / 8.) * P6 * L6 * C6 + (1. / 8.) * P5 * L5 * C5};
  Vector<Real> Nw3{-(1. / 8.) * P4 * L4 * S4 + (1. / 8.) * P6 * L6 * S6,
                   -(1. / 8.) * P5 * L5 * S5 + (1. / 8.) * P4 * L4 * S4,
                   -(1. / 8.) * P6 * L6 * S6 + (1. / 8.) * P5 * L5 * S5};
  Vector<Real> Nx1{3. / (2. * L4) * P4 * C4 - 3. / (2. * L6) * P6 * C6,
                   3. / (2. * L5) * P5 * C5 - 3. / (2. * L4) * P4 * C4,
                   3. / (2. * L6) * P6 * C6 - 3. / (2. * L5) * P5 * C5};
  Vector<Real> Nx2{N1 - (3. / 4.) * P4 * C4 * C4 - (3. / 4.) * P6 * C6 * C6,
                   N2 - (3. / 4.) * P5 * C5 * C5 - (3. / 4.) * P4 * C4 * C4,
                   N3 - (3. / 4.) * P6 * C6 * C6 - (3. / 4.) * P5 * C5 * C5};
  Vector<Real> Nx3{-(3. / 4.) * P4 * C4 * S4 - (3. / 4.) * P6 * C6 * S6,
                   -(3. / 4.) * P5 * C5 * S5 - (3. / 4.) * P4 * C4 * S4,
                   -(3. / 4.) * P6 * C6 * S6 - (3. / 4.) * P5 * C5 * S5};
  Vector<Real> Ny1{3. / (2. * L4) * P4 * S4 - 3. / (2. * L6) * P6 * S6,
                   3. / (2. * L5) * P5 * S5 - 3. / (2. * L4) * P4 * S4,
                   3. / (2. * L6) * P6 * S6 - 3. / (2. * L5) * P5 * S5};
  Vector<Real> Ny2{-(3. / 4.) * P4 * C4 * S4 - (3. / 4.) * P6 * C6 * S6,
                   -(3. / 4.) * P5 * C5 * S5 - (3. / 4.) * P4 * C4 * S4,
                   -(3. / 4.) * P6 * C6 * S6 - (3. / 4.) * P5 * C5 * S5};
  Vector<Real> Ny3{N1 - (3. / 4.) * P4 * S4 * S4 - (3. / 4.) * P6 * S6 * S6,
                   N2 - (3. / 4.) * P5 * S5 * S5 - (3. / 4.) * P4 * S4 * S4,
                   N3 - (3. / 4.) * P6 * S6 * S6 - (3. / 4.) * P5 * S5 * S5};

  // clang-format off
  N = Matrix<Real> {
  //      0    1       2       3       4      5      6       7       8       9     10     11      12      13      14
    {N0(0),    0.,     0.,     0.,     0., N0(1),    0.,     0.,     0.,     0., N0(2),    0.,     0.,     0.,     0.},  // 0
    {   0., N0(0),     0.,     0.,     0.,    0., N0(1),     0.,     0.,     0.,    0., N0(2),     0.,     0.,     0.},  // 1
    {   0.,    0.,  N0(0), Nw2(0), Nw3(0),    0.,    0.,  N0(1), Nw2(1), Nw3(1),    0.,    0.,  N0(2), Nw2(2), Nw3(2)},  // 2
    {   0.,    0., Nx1(0), Nx2(0), Nx3(0),    0.,    0., Nx1(1), Nx2(1), Nx3(1),    0.,    0., Nx1(2), Nx2(2), Nx3(2)},  // 3
    {   0.,    0., Ny1(0), Ny2(0), Ny3(0),    0.,    0., Ny1(1), Ny2(1), Ny3(1),    0.,    0., Ny1(2), Ny2(2), Ny3(2)},  // 4
    {   0,     0.,     0.,     0.,     0.,    0.,    0.,     0.,     0.,     0.,    0.,    0.,     0.,     0.,     0.}}; // 5 ???
  // clang-format on

#endif
}

/* -------------------------------------------------------------------------- */
template <>
inline void
InterpolationElement<_itp_discrete_kirchhoff_triangle_18>::computeDNDS(
    const Vector<Real> & /*natural_coords*/, Matrix<Real> & /*B*/) {

#if 0
  // natural coordinate
  Real xi = natural_coords(0);
  Real eta = natural_coords(1);

  // projected_coord (x1 x2 x3) (y1 y2 y3)

  // donne juste pour pour le patch test 4_5_5 mais donne quelque changement
  // de signe dans la matrice de rotation

  Real x21 = projected_coord(0, 0) - projected_coord(0, 1); // x1-x2
  Real x32 = projected_coord(0, 1) - projected_coord(0, 2);
  Real x13 = projected_coord(0, 2) - projected_coord(0, 0);

  Real y21 = projected_coord(1, 0) - projected_coord(1, 1); // y1-y2
  Real y32 = projected_coord(1, 1) - projected_coord(1, 2);
  Real y13 = projected_coord(1, 2) - projected_coord(1, 0);

  // donne juste pour la matrice de rigidité... mais pas pour le patch test
  // 4_5_5

  /* Real x21=projected_coord(0,1)-projected_coord(0,0);
    Real x32=projected_coord(0,2)-projected_coord(0,1);
    Real x13=projected_coord(0,0)-projected_coord(0,2);

    Real y21=projected_coord(1,1)-projected_coord(1,0);
    Real y32=projected_coord(1,2)-projected_coord(1,1);
    Real y13=projected_coord(1,0)-projected_coord(1,2);*/

  // natural triangle side length
  Real L4 = std::sqrt(x21 * x21 + y21 * y21);
  Real L5 = std::sqrt(x32 * x32 + y32 * y32);
  Real L6 = std::sqrt(x13 * x13 + y13 * y13);

  // sinus and cosinus
  Real C4 = x21 / L4;
  Real C5 = x32 / L5;
  Real C6 = x13 / L6;
  Real S4 = y21 / L4;
  Real S5 = y32 / L5;
  Real S6 = y13 / L6;

  Real dN1xi = -1;
  Real dN2xi = 1;
  Real dN3xi = 0;

  Real dN1eta = -1;
  Real dN2eta = 0;
  Real dN3eta = 1;

  Real dP4xi = 4 - 8. * xi - 4 * eta;
  Real dP5xi = 4 * eta;
  Real dP6xi = -4 * eta;

  Real dP4eta = -4 * xi;
  Real dP5eta = 4 * xi;
  Real dP6eta = 4 - 4 * xi - 8. * eta;

  // N'xi
  auto Np00 = dN1xi;
  auto Np01 = dN2xi;
  auto Np02 = dN3xi;
  //   N'eta
  auto Np10 = dN1eta;
  auto Np11 = dN2eta;
  auto Np12 = dN3eta;

  // Nxi1'xi
  auto Nx1p00 = 3. / (2 * L4) * dP4xi * C4 - 3. / (2. * L6) * dP6xi * C6;
  auto Nx1p01 = 3. / (2 * L5) * dP5xi * C5 - 3. / (2. * L4) * dP4xi * C4;
  auto Nx1p02 = 3. / (2 * L6) * dP6xi * C6 - 3. / (2. * L5) * dP5xi * C5;
  //    Nxi1'eta
  auto Nx1p10 = 3. / (2 * L4) * dP4eta * C4 - 3. / (2. * L6) * dP6eta * C6;
  auto Nx1p11 = 3. / (2 * L5) * dP5eta * C5 - 3. / (2. * L4) * dP4eta * C4;
  auto Nx1p12 = 3. / (2 * L6) * dP6eta * C6 - 3. / (2. * L5) * dP5eta * C5;

  // Nxi2'xi
  auto Nx2p00 = -1 - (3. / 4.) * dP4xi * C4 * C4 - (3. / 4.) * dP6xi * C6 * C6;
  auto Nx2p01 = 1 - (3. / 4.) * dP5xi * C5 * C5 - (3. / 4.) * dP4xi * C4 * C4;
  auto Nx2p02 = -(3. / 4.) * dP6xi * C6 * C6 - (3. / 4.) * dP5xi * C5 * C5;
  //    Nxi2'eta
  auto Nx2p10 =
      -1 - (3. / 4.) * dP4eta * C4 * C4 - (3. / 4.) * dP6eta * C6 * C6;
  auto Nx2p11 = -(3. / 4.) * dP5eta * C5 * C5 - (3. / 4.) * dP4eta * C4 * C4;
  auto Nx2p12 = 1 - (3. / 4.) * dP6eta * C6 * C6 - (3. / 4.) * dP5eta * C5 * C5;

  // Nxi3'xi
  auto Nx3p00 = -(3. / 4.) * dP4xi * C4 * S4 - (3. / 4.) * dP6xi * C6 * S6;
  auto Nx3p01 = -(3. / 4.) * dP5xi * C5 * S5 - (3. / 4.) * dP4xi * C4 * S4;
  auto Nx3p02 = -(3. / 4.) * dP6xi * C6 * S6 - (3. / 4.) * dP5xi * C5 * S5;
  //  Nxi3'eta
  auto Nx3p10 = -(3. / 4.) * dP4eta * C4 * S4 - (3. / 4.) * dP6eta * C6 * S6;
  auto Nx3p11 = -(3. / 4.) * dP5eta * C5 * S5 - (3. / 4.) * dP4eta * C4 * S4;
  auto Nx3p12 = -(3. / 4.) * dP6eta * C6 * S6 - (3. / 4.) * dP5eta * C5 * S5;

  // Nyi1'xi
  auto Ny1p00 = 3 / (2 * L4) * dP4xi * S4 - 3 / (2 * L6) * dP6xi * S6;
  auto Ny1p01 = 3 / (2 * L5) * dP5xi * S5 - 3 / (2 * L4) * dP4xi * S4;
  auto Ny1p02 = 3 / (2 * L6) * dP6xi * S6 - 3 / (2 * L5) * dP5xi * S5;
  //    Nyi1'eta
  auto Ny1p10 = 3 / (2 * L4) * dP4eta * S4 - 3 / (2 * L6) * dP6eta * S6;
  auto Ny1p11 = 3 / (2 * L5) * dP5eta * S5 - 3 / (2 * L4) * dP4eta * S4;
  auto Ny1p12 = 3 / (2 * L6) * dP6eta * S6 - 3 / (2 * L5) * dP5eta * S5;

  // Nyi2'xi
  auto Ny2p00 = -(3. / 4.) * dP4xi * C4 * S4 - (3. / 4.) * dP6xi * C6 * S6;
  auto Ny2p01 = -(3. / 4.) * dP5xi * C5 * S5 - (3. / 4.) * dP4xi * C4 * S4;
  auto Ny2p02 = -(3. / 4.) * dP6xi * C6 * S6 - (3. / 4.) * dP5xi * C5 * S5;
  //  Nyi2'eta
  auto Ny2p10 = -(3. / 4.) * dP4eta * C4 * S4 - (3. / 4.) * dP6eta * C6 * S6;
  auto Ny2p11 = -(3. / 4.) * dP5eta * C5 * S5 - (3. / 4.) * dP4eta * C4 * S4;
  auto Ny2p12 = -(3. / 4.) * dP6eta * C6 * S6 - (3. / 4.) * dP5eta * C5 * S5;

  // Nyi3'xi
  auto Ny3p00 =
      dN1xi - (3. / 4.) * dP4xi * S4 * S4 - (3. / 4.) * dP6xi * S6 * S6;
  auto Ny3p01 =
      dN2xi - (3. / 4.) * dP5xi * S5 * S5 - (3. / 4.) * dP4xi * S4 * S4;
  auto Ny3p02 =
      dN3xi - (3. / 4.) * dP6xi * S6 * S6 - (3. / 4.) * dP5xi * S5 * S5;
  //      Nyi3'eta
  auto Ny3p10 =
      dN1eta - (3. / 4.) * dP4eta * S4 * S4 - (3. / 4.) * dP6eta * S6 * S6;
  auto Ny3p11 =
      dN2eta - (3. / 4.) * dP5eta * S5 * S5 - (3. / 4.) * dP4eta * S4 * S4;
  auto Ny3p12 =
      dN3eta - (3. / 4.) * dP6eta * S6 * S6 - (3. / 4.) * dP5eta * S5 * S5;

  // clang-format off
  //       0    1                2                3                4       5     6                7                8                9    10    11               12               13               14
  B = {{Np00,   0.,              0.,              0.,              0.,   Np01,   0.,              0.,              0.,              0., Np02,   0.,              0.,              0.,              0.},  // 0
       {  0., Np10,              0.,              0.,              0.,     0., Np11,              0.,              0.,              0.,   0., Np12,              0.,              0.,              0.},  // 1
       {Np10, Np00,              0.,              0.,              0.,   Np11, Np01,              0.,              0.,              0., Np12, Np02,              0.,              0.,              0.},  // 2
       {  0.,   0.,          Nx1p00,          Nx2p00,          Nx3p00,     0.,   0.,          Nx1p01,          Nx2p01,          Nx3p01,   0.,   0.,          Nx1p02,          Nx2p02,          Nx3p02},  // 3
       {  0.,   0.,          Ny1p10,          Ny2p10,          Ny3p10,     0.,   0.,          Ny1p11,          Ny2p11,          Ny3p11,   0.,   0.,          Ny1p12,          Ny2p12,          Ny3p12},  // 4
       {  0.,   0., Nx1p10 + Ny1p00, Nx2p10 + Ny2p00, Nx3p10 + Ny3p00,     0.,   0., Nx1p11 + Ny1p01, Nx2p11 + Ny2p01, Nx3p11 + Ny3p01,   0.,   0., Nx1p12 + Ny1p02, Nx2p12 + Ny2p02, Nx3p12 + Ny3p02}}; // 5
  // clang-format on
#endif
}

} // namespace akantu

#endif /* __AKANTU_ELEMENT_CLASS_KIRCHHOFF_SHELL_INLINE_IMPL_CC__ */
