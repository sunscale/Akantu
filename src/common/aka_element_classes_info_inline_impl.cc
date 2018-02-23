/**
 * @file   aka_element_classes_info_inline_impl.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Jun 18 2015
 * @date last modification: Wed Jan 10 2018
 *
 * @brief  Implementation of the streaming fonction for the element classes
 * enums
 *
 * @section LICENSE
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include <unordered_map>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_AKA_ELEMENT_CLASSES_INFO_INLINE_IMPL_CC__
#define __AKANTU_AKA_ELEMENT_CLASSES_INFO_INLINE_IMPL_CC__

AKANTU_ENUM_HASH(ElementType)

#define AKANTU_PP_ELEMTYPE_TO_STR(s, data, elem)                               \
  ({akantu::elem, BOOST_PP_STRINGIZE(elem)})

#define AKANTU_PP_STR_TO_ELEMTYPE(s, data, elem)                               \
  ({BOOST_PP_STRINGIZE(elem), akantu::elem})

namespace aka {
inline std::string to_string(const akantu::ElementType & type) {
  static std::unordered_map<akantu::ElementType, std::string> convert{
      BOOST_PP_SEQ_FOR_EACH_I(
          AKANTU_PP_ENUM, BOOST_PP_SEQ_SIZE(AKANTU_ALL_ELEMENT_TYPE),
          BOOST_PP_SEQ_TRANSFORM(AKANTU_PP_ELEMTYPE_TO_STR, _,
                                 AKANTU_ALL_ELEMENT_TYPE)),
      {akantu::_not_defined, "_not_defined"},
      {akantu::_max_element_type, "_max_element_type"}};
  return convert.at(type);
}
}

namespace akantu {

#define STRINGIFY(type) stream << BOOST_PP_STRINGIZE(type)

/* -------------------------------------------------------------------------- */
//! standard output stream operator for ElementType
inline std::ostream & operator<<(std::ostream & stream,
                                 const ElementType & type) {
  stream << aka::to_string(type);
  return stream;
}

/* -------------------------------------------------------------------------- */

//! standard input stream operator for ElementType
inline std::istream & operator>>(std::istream & stream, ElementType & type) {
  std::string str;
  stream >> str;
  static std::unordered_map<std::string, ElementType> convert{
      BOOST_PP_SEQ_FOR_EACH_I(
          AKANTU_PP_ENUM, BOOST_PP_SEQ_SIZE(AKANTU_ALL_ELEMENT_TYPE),
          BOOST_PP_SEQ_TRANSFORM(AKANTU_PP_STR_TO_ELEMTYPE, _,
                                 AKANTU_ALL_ELEMENT_TYPE))};
  type = convert.at(str);
  return stream;
}

/* -------------------------------------------------------------------------- */
//! standard output stream operator for ElementType
inline std::ostream & operator<<(std::ostream & stream, ElementKind kind) {
  AKANTU_BOOST_ALL_KIND_SWITCH(STRINGIFY);

  return stream;
}

/* -------------------------------------------------------------------------- */
/// standard output stream operator for InterpolationType
inline std::ostream & operator<<(std::ostream & stream,
                                 InterpolationType type) {
  switch (type) {
    BOOST_PP_SEQ_FOR_EACH(AKANTU_BOOST_CASE_MACRO, STRINGIFY,
                          AKANTU_INTERPOLATION_TYPES)
  case _itp_not_defined:
    stream << "_itp_not_defined";
    break;
  }
  return stream;
}

#undef STRINGIFY

} // akantu

#endif /* __AKANTU_AKA_ELEMENT_CLASSES_INFO_INLINE_IMPL_CC__ */
