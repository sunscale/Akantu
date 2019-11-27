/**
 * @file   parser_input_files.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Nov 11 2015
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  implementation of the parser
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

#if defined(__INTEL_COMPILER)
//#pragma warning ( disable : 383 )
#elif defined(__clang__) // test clang to be sure that when we test for gnu it
// is only gnu
#elif (defined(__GNUC__) || defined(__GNUG__))
#define GCC_VERSION                                                            \
  (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#if GCC_VERSION > 40600
#pragma GCC diagnostic push
#endif
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#endif

/* -------------------------------------------------------------------------- */
#include "parser.hh"
#include "parser_grammar_tmpl.hh"
/* -------------------------------------------------------------------------- */
#include "input_file_parser.hh"
/* -------------------------------------------------------------------------- */
#include <fstream>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
void Parser::parse(const std::string & filename) {
  this->clean();
  std::ifstream input(filename.c_str());

  if (!input.good()) {
    AKANTU_EXCEPTION("Could not open file " << filename << "!");
  }

  input.unsetf(std::ios::skipws);

  // wrap istream into iterator
  spirit::istream_iterator fwd_begin(input);
  spirit::istream_iterator fwd_end;

  // wrap forward iterator with position iterator, to record the position
  using pos_iterator_type =
      spirit::classic::position_iterator2<spirit::istream_iterator>;
  pos_iterator_type position_begin(fwd_begin, fwd_end, filename);
  pos_iterator_type position_end;

  // parse
  parser::InputFileGrammar<pos_iterator_type> ag(this);

  bool result = qi::phrase_parse(position_begin, position_end, ag, ag.skipper);

  if (!result || position_begin != position_end) {
    spirit::classic::file_position pos = position_begin.get_position();

    AKANTU_EXCEPTION("Parse error [ "
                     << ag.getErrorMessage() << " ]"
                     << " in file " << filename << " line " << pos.line
                     << " column " << pos.column << std::endl
                     << "'" << position_begin.get_currentline() << "'"
                     << std::endl
                     << std::setw(pos.column) << " "
                     << "^- here");
  }

  try {
    bool permissive = getParameter("permissive_parser", _ppsc_current_scope);

    permissive_parser = permissive;
    AKANTU_DEBUG_INFO("Parser switched permissive mode to "
                      << std::boolalpha << permissive_parser);
  } catch (debug::Exception & e) {
  }

  last_parsed_file = filename;
  input.close();
}

} // namespace akantu

#if defined(__INTEL_COMPILER)
//#pragma warning ( disable : 383 )
#elif defined(__clang__) // test clang to be sure that when we test for gnu it
// is only gnu
#elif defined(__GNUG__)
#if GCC_VERSION > 40600
#pragma GCC diagnostic pop
#else
#pragma GCC diagnostic warning "-Wunused-local-typedefs"
#endif
#endif
