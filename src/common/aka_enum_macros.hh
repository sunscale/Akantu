/**
 * @file   aka_enum_macros.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  Wed Oct 31 2018
 *
 * @brief A Documented file.
 *
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
#include <algorithm>
#include <string>
/* -------------------------------------------------------------------------- */
#ifndef AKANTU_AKA_ENUM_MACROS_HH_
#define AKANTU_AKA_ENUM_MACROS_HH_

#define AKANTU_PP_ENUM(s, data, i, elem)                                       \
  BOOST_PP_TUPLE_REM()                                                         \
  elem BOOST_PP_COMMA_IF(BOOST_PP_NOT_EQUAL(i, BOOST_PP_DEC(data)))

#if (defined(__GNUC__) || defined(__GNUG__))
#define AKA_GCC_VERSION                                                        \
  (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#if AKA_GCC_VERSION < 60000
#define AKANTU_ENUM_HASH(type_name)                                            \
  namespace std {                                                              \
  template <> struct hash<::akantu::type_name> {                               \
    using argument_type = ::akantu::type_name;                                 \
    size_t operator()(const argument_type & e) const noexcept {                \
      auto ue = underlying_type_t<argument_type>(e);                           \
      return uh(ue);                                                           \
    }                                                                          \
                                                                               \
  private:                                                                     \
    const hash<underlying_type_t<argument_type>> uh{};                         \
  };                                                                           \
  }
#else
#define AKANTU_ENUM_HASH(type_name)
#endif // AKA_GCC_VERSION
#endif // GNU

#define AKANTU_PP_CAT(s, data, elem) BOOST_PP_CAT(data, elem)

#define AKANTU_PP_TYPE_TO_STR(s, data, elem)                                   \
  ({BOOST_PP_CAT(data, elem), BOOST_PP_STRINGIZE(elem)})

#define AKANTU_PP_STR_TO_TYPE(s, data, elem)                                   \
  ({BOOST_PP_STRINGIZE(elem), BOOST_PP_CAT(data, elem)})

#define AKANTU_CLASS_ENUM_DECLARE(type_name, list)                             \
  enum class type_name {                                                       \
    BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(AKANTU_PP_CAT, _, list))          \
  };

#define AKANTU_ENUM_OUTPUT_STREAM_(type_name, list, prefix)                    \
  }                                                                            \
  AKANTU_ENUM_HASH(type_name)                                                  \
  namespace std {                                                              \
  inline string to_string(const ::akantu::type_name & type) {                  \
    using namespace akantu;                                                    \
    static unordered_map<::akantu::type_name, string> convert{                 \
        BOOST_PP_SEQ_FOR_EACH_I(                                               \
            AKANTU_PP_ENUM, BOOST_PP_SEQ_SIZE(list),                           \
            BOOST_PP_SEQ_TRANSFORM(AKANTU_PP_TYPE_TO_STR, prefix, list))};     \
    return convert.at(type);                                                   \
  }                                                                            \
  }                                                                            \
  namespace akantu {                                                           \
  inline std::ostream & operator<<(std::ostream & stream,                      \
                                   const type_name & type) {                   \
    stream << std::to_string(type);                                            \
    return stream;                                                             \
  }

#define AKANTU_ENUM_INPUT_STREAM_(type_name, list, prefix)                     \
  inline std::istream & operator>>(std::istream & stream,                      \
                                   type_name & type) { /* NOLINT */            \
    std::string str;                                                           \
    stream >> str; /* NOLINT */                                                \
    static std::unordered_map<std::string, type_name> convert{                 \
        BOOST_PP_SEQ_FOR_EACH_I(                                               \
            AKANTU_PP_ENUM, BOOST_PP_SEQ_SIZE(list),                           \
            BOOST_PP_SEQ_TRANSFORM(AKANTU_PP_STR_TO_TYPE, prefix, list))};     \
    try {                                                                      \
      type = convert.at(str);                                                  \
    } catch (std::out_of_range &) {                                            \
      std::ostringstream values;                                               \
      std::for_each(convert.begin(), convert.end(), [&values](auto && pair) {  \
        static bool first = true;                                              \
        if (not first)                                                         \
          values << ", ";                                                      \
        values << "\"" << pair.first << "\"";                                  \
        first = false;                                                         \
      });                                                                      \
      AKANTU_EXCEPTION("The value " << str << " is not a valid "               \
                                    << BOOST_PP_STRINGIZE(type_name)           \
                                    << " valid values are " << values.str());  \
    }                                                                          \
    return stream;                                                             \
  }

#define AKANTU_CLASS_ENUM_OUTPUT_STREAM(type_name, list)                       \
  AKANTU_ENUM_OUTPUT_STREAM_(type_name, list, type_name::_)

#define AKANTU_ENUM_OUTPUT_STREAM(type_name, list)                             \
  AKANTU_ENUM_OUTPUT_STREAM_(type_name, list, )

#define AKANTU_CLASS_ENUM_INPUT_STREAM(type_name, list)                        \
  AKANTU_ENUM_INPUT_STREAM_(type_name, list, type_name::_)

#define AKANTU_ENUM_INPUT_STREAM(type_name, list)                              \
  AKANTU_ENUM_INPUT_STREAM_(type_name, list, )

#endif /* AKANTU_AKA_ENUM_MACROS_HH_ */
