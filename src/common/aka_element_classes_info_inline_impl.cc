/**
 * @file   aka_element_classes_info_inline_impl.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Jun 18 2015
 * @date last modification: Sun Jul 19 2015
 *
 * @brief  Implementation of the streaming fonction for the element classes enums
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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

namespace akantu {

#define STRINGIFY(type)				\
  stream << BOOST_PP_STRINGIZE(type)

/* -------------------------------------------------------------------------- */
//! standard output stream operator for ElementType
inline std::ostream & operator <<(std::ostream & stream, ElementType type) {
  switch(type) {
    BOOST_PP_SEQ_FOR_EACH(AKANTU_BOOST_CASE_MACRO, \
                          STRINGIFY,               \
                          AKANTU_ALL_ELEMENT_TYPE)
    case _not_defined:       stream << "_not_defined"; break;
    case _max_element_type:  stream << "_max_element_type"; break;
  }

  return stream;
}

/* -------------------------------------------------------------------------- */

//! standard input stream operator for ElementType
inline std::istream & operator >>(std::istream & stream, ElementType & type) {
#define IF_SEQUENCE(_type)				\
  else if (tmp == BOOST_PP_STRINGIZE(_type)) type = _type;

  std::string tmp;
  stream >> tmp;

  if (1 == 2) {}
  BOOST_PP_SEQ_FOR_EACH(AKANTU_BOOST_LIST_MACRO,     \
			IF_SEQUENCE,		     \
			AKANTU_ALL_ELEMENT_TYPE)
  else AKANTU_EXCEPTION("unknown element type: '" << tmp << "'");

#undef IF_SEQUENCE
  return stream;
}

/* -------------------------------------------------------------------------- */
//! standard output stream operator for ElementType
inline std::ostream & operator <<(std::ostream & stream, ElementKind kind ) {
  AKANTU_BOOST_ALL_KIND_SWITCH(STRINGIFY);

  return stream;
}

/* -------------------------------------------------------------------------- */
/// standard output stream operator for InterpolationType
inline std::ostream & operator <<(std::ostream & stream, InterpolationType type)
{
  switch(type) {
    BOOST_PP_SEQ_FOR_EACH(AKANTU_BOOST_CASE_MACRO, \
                          STRINGIFY,               \
                          AKANTU_INTERPOLATION_TYPES)
    case _itp_not_defined             : stream << "_itp_not_defined"            ; break;
    }
  return stream;
}

#undef STRINGIFY

} // akantu
