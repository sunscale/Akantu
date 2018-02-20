/**
 * @file   aka_common_inline_impl.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Thu Oct 15 2015
 *
 * @brief  inline implementations of common akantu type descriptions
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
 * All common things to be included in the projects files
 *
 */

/* -------------------------------------------------------------------------- */
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <iostream>
#include "aka_common.hh"

/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_AKA_COMMON_INLINE_IMPL_CC__
#define __AKANTU_AKA_COMMON_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
/// standard output stream operator for GhostType
inline std::ostream & operator<<(std::ostream & stream, GhostType type) {
  switch (type) {
  case _not_ghost:
    stream << "not_ghost";
    break;
  case _ghost:
    stream << "ghost";
    break;
  case _casper:
    stream << "Casper the friendly ghost";
    break;
  }
  return stream;
}

/* -------------------------------------------------------------------------- */
/// standard output stream operator for TimeStepSolverType
inline std::ostream & operator<<(std::ostream & stream,
                                 const TimeStepSolverType & type) {
  switch (type) {
  case _tsst_static:
    stream << "static";
    break;
  case _tsst_dynamic:
    stream << "dynamic";
    break;
  case _tsst_dynamic_lumped:
    stream << "dynamic_lumped";
    break;
  case _tsst_not_defined:
    stream << "not defined time step solver";
    break;
  }
  return stream;
}

/* -------------------------------------------------------------------------- */
/// standard input stream operator for TimeStepSolverType
inline std::istream & operator>>(std::istream & stream,
                                 TimeStepSolverType & type) {
  std::string str;
  stream >> str;
  if (str == "static")
    type = _tsst_static;
  else if (str == "dynamic")
    type = _tsst_dynamic;
  else if (str == "dynamic_lumped")
    type = _tsst_dynamic_lumped;
  else {
    AKANTU_ERROR("The type "
                       << str << " is not a recognized TimeStepSolverType");

    stream.setstate(std::ios::failbit);
  }

  return stream;
}

/* -------------------------------------------------------------------------- */
/// standard output stream operator for NonLinearSolverType
inline std::ostream & operator<<(std::ostream & stream,
                                 const NonLinearSolverType & type) {
  switch (type) {
  case _nls_linear:
    stream << "linear";
    break;
  case _nls_newton_raphson:
    stream << "newton_raphson";
    break;
  case _nls_newton_raphson_modified:
    stream << "newton_raphson_modified";
    break;
  case _nls_lumped:
    stream << "lumped";
    break;
  case _nls_auto:
    stream << "auto";
    break;
  }
  return stream;
}

/* -------------------------------------------------------------------------- */
/// standard input stream operator for NonLinearSolverType
inline std::istream & operator>>(std::istream & stream,
                                 NonLinearSolverType & type) {
  std::string str;
  stream >> str;
  if (str == "linear")
    type = _nls_linear;
  else if (str == "newton_raphson")
    type = _nls_newton_raphson;
  else if (str == "newton_raphson_modified")
    type = _nls_newton_raphson_modified;
  else if (str == "lumped")
    type = _nls_lumped;
  else if (str == "auto")
    type = _nls_auto;
  else
    type = _nls_auto;

  return stream;
}

/* -------------------------------------------------------------------------- */
/// standard output stream operator for IntegrationSchemeType
inline std::ostream & operator<<(std::ostream & stream,
                                 const IntegrationSchemeType & type) {
  switch (type) {
  case _ist_pseudo_time:
    stream << "pseudo_time";
    break;
  case _ist_forward_euler:
    stream << "forward_euler";
    break;
  case _ist_trapezoidal_rule_1:
    stream << "trapezoidal_rule_1";
    break;
  case _ist_backward_euler:
    stream << "backward_euler";
    break;
  case _ist_central_difference:
    stream << "central_difference";
    break;
  case _ist_fox_goodwin:
    stream << "fox_goodwin";
    break;
  case _ist_trapezoidal_rule_2:
    stream << "trapezoidal_rule_2";
    break;
  case _ist_linear_acceleration:
    stream << "linear_acceleration";
    break;
  case _ist_newmark_beta:
    stream << "newmark_beta";
    break;
  case _ist_generalized_trapezoidal:
    stream << "generalized_trapezoidal";
    break;
  }
  return stream;
}

/* -------------------------------------------------------------------------- */
/// standard input stream operator for IntegrationSchemeType
inline std::istream & operator>>(std::istream & stream,
                                 IntegrationSchemeType & type) {
  std::string str;
  stream >> str;

  if (str == "pseudo_time")
    type = _ist_pseudo_time;
  else if (str == "forward_euler")
    type = _ist_forward_euler;
  else if (str == "trapezoidal_rule_1")
    type = _ist_trapezoidal_rule_1;
  else if (str == "backward_euler")
    type = _ist_backward_euler;
  else if (str == "central_difference")
    type = _ist_central_difference;
  else if (str == "fox_goodwin")
    type = _ist_fox_goodwin;
  else if (str == "trapezoidal_rule_2")
    type = _ist_trapezoidal_rule_2;
  else if (str == "linear_acceleration")
    type = _ist_linear_acceleration;
  else if (str == "newmark_beta")
    type = _ist_newmark_beta;
  else if (str == "generalized_trapezoidal")
    type = _ist_generalized_trapezoidal;
  else {
    AKANTU_ERROR("The type "
                       << str << " is not a recognized IntegrationSchemeType");
    stream.setstate(std::ios::failbit);
  }
  return stream;
}

/* -------------------------------------------------------------------------- */
/// standard output stream operator for SynchronizationTag
inline std::ostream & operator<<(std::ostream & stream,
                                 SynchronizationTag type) {
  switch (type) {
  case _gst_whatever:
    stream << "_gst_whatever";
    break;
  case _gst_update:
    stream << "_gst_update";
    break;
  case _gst_size:
    stream << "_gst_size";
    break;
  case _gst_smm_mass:
    stream << "_gst_smm_mass";
    break;
  case _gst_smm_for_gradu:
    stream << "_gst_smm_for_gradu";
    break;
  case _gst_smm_boundary:
    stream << "_gst_smm_boundary";
    break;
  case _gst_smm_uv:
    stream << "_gst_smm_uv";
    break;
  case _gst_smm_res:
    stream << "_gst_smm_res";
    break;
  case _gst_smm_init_mat:
    stream << "_gst_smm_init_mat";
    break;
  case _gst_smm_stress:
    stream << "_gst_smm_stress";
    break;
  case _gst_smmc_facets:
    stream << "_gst_smmc_facets";
    break;
  case _gst_smmc_facets_conn:
    stream << "_gst_smmc_facets_conn";
    break;
  case _gst_smmc_facets_stress:
    stream << "_gst_smmc_facets_stress";
    break;
  case _gst_smmc_damage:
    stream << "_gst_smmc_damage";
    break;
  case _gst_giu_global_conn:
    stream << "_gst_giu_global_conn";
    break;
  case _gst_ce_groups:
    stream << "_gst_ce_groups";
    break;
  case _gst_gm_clusters:
    stream << "_gst_gm_clusters";
    break;
  case _gst_htm_capacity:
    stream << "_gst_htm_capacity";
    break;
  case _gst_htm_temperature:
    stream << "_gst_htm_temperature";
    break;
  case _gst_htm_gradient_temperature:
    stream << "_gst_htm_gradient_temperature";
    break;
  case _gst_htm_phi:
    stream << "_gst_htm_phi";
    break;
  case _gst_htm_gradient_phi:
    stream << "_gst_htm_gradient_phi";
    break;
  case _gst_mnl_for_average:
    stream << "_gst_mnl_for_average";
    break;
  case _gst_mnl_weight:
    stream << "_gst_mnl_weight";
    break;
  case _gst_nh_criterion:
    stream << "_gst_nh_criterion";
    break;
  case _gst_test:
    stream << "_gst_test";
    break;
  case _gst_user_1:
    stream << "_gst_user_1";
    break;
  case _gst_user_2:
    stream << "_gst_user_2";
    break;
  case _gst_material_id:
    stream << "_gst_material_id";
    break;
  case _gst_for_dump:
    stream << "_gst_for_dump";
    break;
  case _gst_cf_nodal:
    stream << "_gst_cf_nodal";
    break;
  case _gst_cf_incr:
    stream << "_gst_cf_incr";
    break;
  case _gst_solver_solution:
    stream << "_gst_solver_solution";
    break;
  }
  return stream;
}

/* -------------------------------------------------------------------------- */
/// standard output stream operator for SolveConvergenceCriteria
inline std::ostream & operator<<(std::ostream & stream,
                                 const SolveConvergenceCriteria & criteria) {
  switch (criteria) {
  case _scc_residual:
    stream << "_scc_residual";
    break;
  case _scc_solution:
    stream << "_scc_solution";
    break;
  case _scc_residual_mass_wgh:
    stream << "_scc_residual_mass_wgh";
    break;
  }
  return stream;
}

inline std::istream & operator>>(std::istream & stream,
                                 SolveConvergenceCriteria & criteria) {
  std::string str;
  stream >> str;
  if (str == "residual")
    criteria = _scc_residual;
  else if (str == "solution")
    criteria = _scc_solution;
  else if (str == "residual_mass_wgh")
    criteria = _scc_residual_mass_wgh;
  else {
    stream.setstate(std::ios::failbit);
  }
  return stream;
}

/* -------------------------------------------------------------------------- */
inline std::string to_lower(const std::string & str) {
  std::string lstr = str;
  std::transform(lstr.begin(), lstr.end(), lstr.begin(), (int (*)(int))tolower);
  return lstr;
}

/* -------------------------------------------------------------------------- */
inline std::string trim(const std::string & to_trim) {
  std::string trimed = to_trim;
  // left trim
  trimed.erase(trimed.begin(),
               std::find_if(trimed.begin(), trimed.end(),
                            std::not1(std::ptr_fun<int, int>(isspace))));
  // right trim
  trimed.erase(std::find_if(trimed.rbegin(), trimed.rend(),
                            std::not1(std::ptr_fun<int, int>(isspace)))
                   .base(),
               trimed.end());
  return trimed;
}

} // akantu

#include <cmath>

namespace akantu {

/* -------------------------------------------------------------------------- */
template <typename T> std::string printMemorySize(UInt size) {
  Real real_size = size * sizeof(T);

  UInt mult = 0;
  if (real_size != 0)
    mult = (std::log(real_size) / std::log(2)) / 10;

  std::stringstream sstr;

  real_size /= Real(1 << (10 * mult));
  sstr << std::setprecision(2) << std::fixed << real_size;

  std::string size_prefix;
  switch (mult) {
  case 0:
    sstr << "";
    break;
  case 1:
    sstr << "Ki";
    break;
  case 2:
    sstr << "Mi";
    break;
  case 3:
    sstr << "Gi";
    break; // I started on this type of machines
           // (32bit computers) (Nicolas)
  case 4:
    sstr << "Ti";
    break;
  case 5:
    sstr << "Pi";
    break;
  case 6:
    sstr << "Ei";
    break; // theoritical limit of RAM of the current
           // computers in 2014 (64bit computers) (Nicolas)
  case 7:
    sstr << "Zi";
    break;
  case 8:
    sstr << "Yi";
    break;
  default:
    AKANTU_ERROR(
        "The programmer in 2014 didn't thought so far (even wikipedia does not "
        "go further)."
        << " You have at least 1024 times more than a yobibit of RAM!!!"
        << " Just add the prefix corresponding in this switch case.");
  }

  sstr << "Byte";

  return sstr.str();
}

} // akantu

#endif /* __AKANTU_AKA_COMMON_INLINE_IMPL_CC__ */
