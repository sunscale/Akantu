/**
 * @file   stress_based_weight_function_inline_impl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Cyprien Wolff <cyprien.wolff@epfl.ch>
 *
 * @date creation: Fri Apr 13 2012
 * @date last modification: Wed Feb 03 2016
 *
 * @brief  Implementation of inline function of remove damaged with
 * damage rate weight function
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */
/* -------------------------------------------------------------------------- */
#include "stress_based_weight_function.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
/* -------------------------------------------------------------------------- */
inline void StressBasedWeightFunction::updateInternals() {
  // updatePrincipalStress(_not_ghost);
  // updatePrincipalStress(_ghost);
}

/* -------------------------------------------------------------------------- */
// inline void StressBasedWeightFunction::selectType(ElementType type1,
// 						  GhostType ghost_type1,
// 						  ElementType type2,
// 						  GhostType ghost_type2) {
//   selected_stress_diag = &stress_diag(type2, ghost_type2);
//   selected_stress_base = &stress_base(type2, ghost_type2);

//   selected_characteristic_size = &characteristic_size(type1, ghost_type1);
// }

/* -------------------------------------------------------------------------- */
inline Real StressBasedWeightFunction::computeRhoSquare(
    __attribute__((unused)) Real r, __attribute__((unused)) Vector<Real> & eigs,
    __attribute__((unused)) Matrix<Real> & eigenvects,
    __attribute__((unused)) Vector<Real> & x_s) {
  //   if (spatial_dimension == 1)
  //     return eigs[0];
  //   else if (spatial_dimension == 2) {
  //     Vector<Real> u1(eigenvects.storage(), 2);
  //     Real cos_t = x_s.dot(u1) / (x_s.norm() * u1.norm());

  //     Real cos_t_2;
  //     Real sin_t_2;

  //     Real sigma1_2 = eigs[0]*eigs[0];
  //     Real sigma2_2 = eigs[1]*eigs[1];

  // #ifdef __trick__
  //     Real zero = std::numeric_limits<Real>::epsilon();
  //     if(std::abs(cos_t) < zero) {
  //       cos_t_2 = 0;
  //       sin_t_2 = 1;
  //     } else {
  //       cos_t_2 = cos_t * cos_t;
  //       sin_t_2 = (1 - cos_t_2);
  //     }

  //     Real rhop1 = std::max(0., cos_t_2 / sigma1_2);
  //     Real rhop2 = std::max(0., sin_t_2 / sigma2_2);
  // #else
  //     cos_t_2 = cos_t * cos_t;
  //     sin_t_2 = (1 - cos_t_2);

  //     Real rhop1 = cos_t_2 / sigma1_2;
  //     Real rhop2 = sin_t_2 / sigma2_2;
  // #endif

  //     return 1./ (rhop1 + rhop2);
  //   } else if (spatial_dimension == 3) {
  //     Vector<Real> u1(eigenvects.storage() + 0*3, 3);
  //     //Vector<Real> u2(eigenvects.storage() + 1*3, 3);
  //     Vector<Real> u3(eigenvects.storage() + 2*3, 3);

  //     Real zero = std::numeric_limits<Real>::epsilon();

  //     Vector<Real> tmp(3);
  //     tmp.crossProduct(x_s, u3);

  //     Vector<Real> u3_C_x_s_C_u3(3);
  //     u3_C_x_s_C_u3.crossProduct(u3, tmp);

  //     Real norm_u3_C_x_s_C_u3 = u3_C_x_s_C_u3.norm();
  //     Real cos_t = 0.;
  //     if(std::abs(norm_u3_C_x_s_C_u3) > zero) {
  //       Real inv_norm_u3_C_x_s_C_u3 = 1. / norm_u3_C_x_s_C_u3;
  //       cos_t = u1.dot(u3_C_x_s_C_u3) * inv_norm_u3_C_x_s_C_u3;
  //     }

  //     Real cos_p = u3.dot(x_s) / r;

  //     Real cos_t_2;
  //     Real sin_t_2;
  //     Real cos_p_2;
  //     Real sin_p_2;

  //     Real sigma1_2 = eigs[0]*eigs[0];
  //     Real sigma2_2 = eigs[1]*eigs[1];
  //     Real sigma3_2 = eigs[2]*eigs[2];

  // #ifdef __trick__
  //     if(std::abs(cos_t) < zero) {
  //       cos_t_2 = 0;
  //       sin_t_2 = 1;
  //     } else {
  //       cos_t_2 = cos_t * cos_t;
  //       sin_t_2 = (1 - cos_t_2);
  //     }

  //     if(std::abs(cos_p) < zero) {
  //       cos_p_2 = 0;
  //       sin_p_2 = 1;
  //     } else {
  //       cos_p_2 = cos_p * cos_p;
  //       sin_p_2 = (1 - cos_p_2);
  //     }

  //     Real rhop1 = std::max(0., sin_p_2 * cos_t_2 / sigma1_2);
  //     Real rhop2 = std::max(0., sin_p_2 * sin_t_2 / sigma2_2);
  //     Real rhop3 = std::max(0., cos_p_2 / sigma3_2);
  // #else
  //     cos_t_2 = cos_t * cos_t;
  //     sin_t_2 = (1 - cos_t_2);

  //     cos_p_2 = cos_p * cos_p;
  //     sin_p_2 = (1 - cos_p_2);

  //     Real rhop1 = sin_p_2 * cos_t_2 / sigma1_2;
  //     Real rhop2 = sin_p_2 * sin_t_2 / sigma2_2;
  //     Real rhop3 = cos_p_2 / sigma3_2;
  // #endif

  //     return 1./ (rhop1 + rhop2 + rhop3);
  //   }
  return 0.;
}

/* -------------------------------------------------------------------------- */
inline Real StressBasedWeightFunction::operator()(__attribute__((unused))
                                                  Real r,
                                                  __attribute__((unused))
                                                  const IntegrationPoint & q1,
                                                  __attribute__((unused))
                                                  const IntegrationPoint & q2) {
  // Real zero = std::numeric_limits<Real>::epsilon();

  // if(r < zero) return 1.; // means x and s are the same points

  // const Vector<Real> & x = q1.getPosition();
  // const Vector<Real> & s = q2.getPosition();

  // Vector<Real> eigs =
  //   selected_stress_diag->begin(spatial_dimension)[q2.global_num];

  // Matrix<Real> eigenvects =
  //   selected_stress_base->begin(spatial_dimension,
  //   spatial_dimension)[q2.global_num];

  // Real min_rho_lc = selected_characteristic_size->begin()[q1.global_num];

  // Vector<Real> x_s(spatial_dimension);
  // x_s  = x;
  // x_s -= s;

  // Real rho_2 = computeRhoSquare(r, eigs, eigenvects, x_s);

  // Real rho_lc_2 = std::max(this->R2 * rho_2, min_rho_lc*min_rho_lc);

  // // Real w = std::max(0., 1. - r*r / rho_lc_2);
  // // w = w*w;
  // Real w = exp(- 2*2*r*r / rho_lc_2);
  // return w;
  return 0.;
}
} // namespace akantu
