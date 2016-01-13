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
RandomParameter<Real> Parser::parseRandomParameter(const std::string & value, const ParserSection & section) {
  using boost::spirit::ascii::space_type;
  parser::RandomGeneratorGrammar<std::string::const_iterator, space_type> grammar(section);
  grammar.name("random_grammar");
  parser::ParsableRandomGenerator rg = Parser::parseType<parser::ParsableRandomGenerator>(value, grammar);
  Vector<Real> params = rg.parameters;
  switch(rg.type) {
    case _rdt_not_defined: return RandomParameter<Real>(rg.base);
    case _rdt_uniform:     return RandomParameter<Real>(rg.base, UniformDistribution<Real>(params(0), params(1)));
    case _rdt_weibull:     return RandomParameter<Real>(rg.base, WeibullDistribution<Real>(params(0), params(1)));
    default:
      AKANTU_EXCEPTION("This is an unknown random distribution in the parser");
  }
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__
