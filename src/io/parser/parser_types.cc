/**
 * @file   parser.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Nov 13 2013
 * @date last modification: Fri Sep 05 2014
 *
 * @brief  implementation of the parser
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
#include "parser.hh"
#include "parser_grammar_tmpl.hh"
/* -------------------------------------------------------------------------- */
#include "algebraic_parser.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
Vector<Real> Parser::parseVector(const std::string & value, const ParserSection & section) {
  using boost::spirit::ascii::space_type;
  parser::VectorGrammar<std::string::const_iterator, space_type> grammar(section);
  grammar.name("vector_grammar");
  return Parser::parseType<parser::parsable_vector>(value, grammar);
}

/* -------------------------------------------------------------------------- */
Matrix<Real> Parser::parseMatrix(const std::string & value, const ParserSection & section) {
  using boost::spirit::ascii::space_type;
  parser::MatrixGrammar<std::string::const_iterator, space_type> grammar(section);
  grammar.name("matrix_grammar");
  return Parser::parseType<parser::parsable_matrix>(value, grammar);
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__