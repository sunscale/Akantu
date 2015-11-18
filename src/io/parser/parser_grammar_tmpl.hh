/**
 * @file   parsable_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Nov 13 2013
 * @date last modification: Tue Jun 24 2014
 *
 * @brief  implementation of the templated part of ParsableParam Parsable and
 * ParsableParamTyped
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
//#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/support_multi_pass.hpp>
#include <boost/spirit/include/classic_position_iterator.hpp>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_PARSER_GRAMMAR_TMPL_HH
#define AKANTU_PARSER_GRAMMAR_TMPL_HH

__BEGIN_AKANTU__

namespace qi = boost::spirit::qi;


/* -------------------------------------------------------------------------- */
template <class T, class Grammar>
T Parser::parseType(const std::string & value, Grammar & grammar) {
  using boost::spirit::ascii::space;

  std::string::const_iterator b = value.begin();
  std::string::const_iterator e = value.end();

  T resultat = T();
  bool res = false;
  try {
    res = qi::phrase_parse(b, e, grammar, space, resultat);
  } catch(debug::Exception & ex) {
    AKANTU_EXCEPTION("Could not parse '" << value << "' as a "
                     << debug::demangle(typeid(T).name()) << ", an unknown error append '"
                     << ex.what());
  }

  if(!res || (b != e)) {
    AKANTU_EXCEPTION("Could not parse '" << value << "' as a "
                     << debug::demangle(typeid(T).name()) << ", an unknown error append '"
                     << std::string(value.begin(), b) << "<HERE>" << std::string(b, e) << "'");
  }
  return resultat;
}

namespace parser {
  template<class Iterator>
  struct Skipper {
    typedef qi::rule<Iterator, void()> type;
  };
}

__END_AKANTU__

#endif //AKANTU_PARSER_GRAMMAR_TMPL_HH