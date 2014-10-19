/**
 * @file   aka_common_inline_impl.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Dec 01 2011
 * @date last modification: Wed Jul 23 2014
 *
 * @brief  inline implementations of common akantu type descriptions
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
 * All common things to be included in the projects files
 *
 */

#include <algorithm>
#include <iomanip>

#include <cctype>

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
//! standard output stream operator for ElementType
inline std::ostream & operator <<(std::ostream & stream, ElementType type)
{
#define STRINGIFY(type)				\
  stream << BOOST_PP_STRINGIZE(type)

  switch(type) {
    BOOST_PP_SEQ_FOR_EACH(AKANTU_BOOST_CASE_MACRO, \
                          STRINGIFY,               \
                          AKANTU_ALL_ELEMENT_TYPE)
    case _not_defined:       stream << "_not_defined"; break;
    case _max_element_type:  stream << "_max_element_type"; break;
  }



#undef STRINGIFY
  return stream;
}

/* -------------------------------------------------------------------------- */
//! standard output stream operator for ElementType
inline std::ostream & operator <<(std::ostream & stream, ElementKind kind )
{
#define STRINGIFY(kind)				\
  stream << BOOST_PP_STRINGIZE(kind)

  AKANTU_BOOST_ALL_KIND_SWITCH(STRINGIFY);
#undef STRINGIFY
  return stream;
}

/* -------------------------------------------------------------------------- */
/// standard output stream operator for InterpolationType
inline std::ostream & operator <<(std::ostream & stream, InterpolationType type)
{
  switch(type)
    {
    case _itp_lagrange_point_1        : stream << "_itp_lagrange_point_1"       ; break;
    case _itp_lagrange_segment_2      : stream << "_itp_lagrange_segment_2"     ; break;
    case _itp_lagrange_segment_3      : stream << "_itp_lagrange_segment_3"     ; break;
    case _itp_lagrange_triangle_3     : stream << "_itp_lagrange_triangle_3"    ; break;
    case _itp_lagrange_triangle_6     : stream << "_itp_lagrange_triangle_6"    ; break;
    case _itp_lagrange_quadrangle_4   : stream << "_itp_lagrange_quadrangle_4"  ; break;
    case _itp_serendip_quadrangle_8   : stream << "_itp_serendip_quadrangle_8"  ; break;
    case _itp_lagrange_tetrahedron_4  : stream << "_itp_lagrange_tetrahedron_4" ; break;
    case _itp_lagrange_tetrahedron_10 : stream << "_itp_lagrange_tetrahedron_10"; break;
    case _itp_lagrange_hexahedron_8   : stream << "_itp_lagrange_hexahedron_8"  ; break;
    case _itp_lagrange_pentahedron_6  : stream << "_itp_lagrange_pentahedron_6" ; break;
#if defined(AKANTU_STRUCTURAL_MECHANICS)
    case _itp_bernoulli_beam          : stream << "_itp_bernoulli_beam"         ; break;
    case _itp_kirchhoff_shell         : stream << "_itp_kirchhoff_shell"        ; break;
#endif
    case _itp_not_defined             : stream << "_itp_not_defined"            ; break;
    }
  return stream;
}

/* -------------------------------------------------------------------------- */
/// standard output stream operator for GhostType
inline std::ostream & operator <<(std::ostream & stream, GhostType type)
{
  switch(type)
    {
    case _not_ghost : stream << "not_ghost"; break;
    case _ghost     : stream << "ghost"    ; break;
    case _casper    : stream << "Casper the friendly ghost"; break;
    }
  return stream;
}

/* -------------------------------------------------------------------------- */
/// standard output stream operator for SynchronizationTag
inline std::ostream & operator <<(std::ostream & stream, SynchronizationTag type)
{
  switch(type)
    {
    case _gst_smm_mass                 : stream << "_gst_smm_mass"                ; break;
    case _gst_smm_for_gradu            : stream << "_gst_smm_for_gradu"           ; break;
    case _gst_smm_boundary             : stream << "_gst_smm_boundary"            ; break;
    case _gst_smm_uv                   : stream << "_gst_smm_uv"                  ; break;
    case _gst_smm_res                  : stream << "_gst_smm_res"                 ; break;
    case _gst_smm_init_mat             : stream << "_gst_smm_init_mat"            ; break;
    case _gst_smm_stress               : stream << "_gst_smm_stress"              ; break;
    case _gst_smmc_facets              : stream << "_gst_smmc_facets"             ; break;
    case _gst_smmc_facets_conn         : stream << "_gst_smmc_facets_conn"        ; break;
    case _gst_smmc_facets_stress       : stream << "_gst_smmc_facets_stress"      ; break;
    case _gst_smmc_damage              : stream << "_gst_smmc_damage"             ; break;
    case _gst_ce_inserter              : stream << "_gst_ce_inserter"             ; break;
    case _gst_gm_clusters              : stream << "_gst_gm_clusters"             ; break;
    case _gst_htm_capacity             : stream << "_gst_htm_capacity"            ; break;
    case _gst_htm_temperature          : stream << "_gst_htm_temperature"         ; break;
    case _gst_htm_gradient_temperature : stream << "_gst_htm_gradient_temperature"; break;
    case _gst_htm_phi                  : stream << "_gst_htm_phi"                 ; break;
    case _gst_htm_gradient_phi         : stream << "_gst_htm_gradient_phi"        ; break;
    case _gst_mnl_for_average          : stream << "_gst_mnl_for_average"         ; break;
    case _gst_mnl_weight               : stream << "_gst_mnl_weight"              ; break;
    case _gst_test                     : stream << "_gst_test"                    ; break;
    case _gst_material_id              : stream << "_gst_material_id"             ; break;
    case _gst_for_dump                 : stream << "_gst_for_dump"                ; break;
    case _gst_cf_nodal                 : stream << "_gst_cf_nodal"                ; break;
    case _gst_cf_incr                  : stream << "_gst_cf_incr"                 ; break;
    }
  return stream;
}

/* -------------------------------------------------------------------------- */
/// standard output stream operator for SolveConvergenceCriteria
inline std::ostream & operator <<(std::ostream & stream, SolveConvergenceCriteria criteria)
{
  switch(criteria) {
  case _scc_residual : stream << "_scc_residual" ; break;
  case _scc_increment: stream << "_scc_increment"; break;
  case _scc_residual_mass_wgh: stream << "_scc_residual_mass_wgh"; break;
  }
  return stream;
}


/* -------------------------------------------------------------------------- */
inline std::string to_lower(const std::string & str) {
  std::string lstr = str;
  std::transform(lstr.begin(),
                 lstr.end(),
                 lstr.begin(),
                 (int(*)(int))tolower);
  return lstr;
}

/* -------------------------------------------------------------------------- */
inline std::string trim(const std::string & to_trim) {
  std::string trimed = to_trim;
  //left trim
  trimed.erase(trimed.begin(),
               std::find_if(trimed.begin(),
                            trimed.end(),
                            std::not1(std::ptr_fun<int, int>(isspace))));
  // right trim
  trimed.erase(std::find_if(trimed.rbegin(),
                            trimed.rend(),
                            std::not1(std::ptr_fun<int, int>(isspace))).base(),
               trimed.end());
  return trimed;
}

__END_AKANTU__

#include <cmath>

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<typename T>
std::string printMemorySize(UInt size) {
  Real real_size = size * sizeof(T);

  UInt mult = (std::log(real_size) / std::log(2)) / 10;

  std::stringstream sstr;

  real_size /= Real(1 << (10 * mult));
  sstr << std::setprecision(2) << std::fixed << real_size;

  std::string size_prefix;
  switch(mult) {
  case 0: sstr << "";   break;
  case 1: sstr << "Ki"; break;
  case 2: sstr << "Mi"; break;
  case 3: sstr << "Gi"; break; // I started on this type of machines
                               // (32bit computers) (Nicolas)
  case 4: sstr << "Ti"; break;
  case 5: sstr << "Pi"; break;
  case 6: sstr << "Ei"; break; // theoritical limit of RAM of the current
                               // computers in 2014 (64bit computers) (Nicolas)
  case 7: sstr << "Zi"; break;
  case 8: sstr << "Yi"; break;
  default:
    AKANTU_DEBUG_ERROR("The programmer in 2014 didn't thought so far (even wikipedia does not go further)."
                       << " You have at least 1024 times more than a yobibit of RAM!!!"
                       << " Just add the prefix corresponding in this switch case.");
  }

  sstr << "Byte";

  return sstr.str();
}



__END_AKANTU__
