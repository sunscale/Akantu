/**
 * @file   aka_element_classes_info_inline_impl.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Tue May 19 11:48:58 2015
 *
 * @brief  Implementation of the streaming fonction for the element classes enums
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
    case _itp_serendip_hexahedron_20  : stream << "_itp_serendip_hexahedron_20" ; break;
    case _itp_lagrange_pentahedron_15 : stream << "_itp_lagrange_pentahedron_15"; break;
#if defined(AKANTU_STRUCTURAL_MECHANICS)
    case _itp_bernoulli_beam          : stream << "_itp_bernoulli_beam"         ; break;
    case _itp_kirchhoff_shell         : stream << "_itp_kirchhoff_shell"        ; break;
#endif
#if defined(AKANTU_IGFEM)
    case _itp_igfem_segment_3         : stream << "_itp_igfem_segment_3"        ; break;
    case _itp_igfem_triangle_4        : stream << "_itp_igfem_triangle_4"       ; break;
    case _itp_igfem_triangle_5        : stream << "_itp_igfem_triangle_5"       ; break;
#endif
    case _itp_not_defined             : stream << "_itp_not_defined"            ; break;
    }
  return stream;
}

__END_AKANTU__
