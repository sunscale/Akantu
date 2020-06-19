/**
 * @file   parser_random.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Sun Jul 09 2017
 *
 * @brief  implementation of the parser
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#if !defined(DOXYGEN)
#include "parser_grammar_tmpl.hh"
/* -------------------------------------------------------------------------- */
#include "algebraic_parser.hh"
#endif
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
RandomParameter<Real>
Parser::parseRandomParameter(const std::string & value,
                             const ParserSection & section) {
#if !defined(DOXYGEN)
  using boost::spirit::ascii::space_type;
  parser::RandomGeneratorGrammar<std::string::const_iterator, space_type>
      grammar(section);
  grammar.name("random_grammar");
  parser::ParsableRandomGenerator rg =
      Parser::parseType<parser::ParsableRandomGenerator>(value, grammar);
  Vector<Real> params = rg.parameters;

  switch (rg.type) {
  case _rdt_not_defined:
    return RandomParameter<Real>(rg.base,
                                 std::uniform_real_distribution<Real>(0, 0));
  case _rdt_uniform:
    return RandomParameter<Real>(
        rg.base, std::uniform_real_distribution<Real>(params(0), params(1)));
  case _rdt_exponential:
    return RandomParameter<Real>(
        rg.base, std::exponential_distribution<Real>(params(0)));
  case _rdt_gamma:
    return RandomParameter<Real>(
        rg.base, std::gamma_distribution<Real>(params(0), params(1)));
  case _rdt_weibull:
    return RandomParameter<Real>(
        rg.base, std::weibull_distribution<Real>(params(1), params(0)));
  case _rdt_extreme_value:
    return RandomParameter<Real>(
        rg.base, std::extreme_value_distribution<Real>(params(0), params(1)));
  case _rdt_normal:
    return RandomParameter<Real>(
        rg.base, std::normal_distribution<Real>(params(0), params(1)));
  case _rdt_lognormal:
    return RandomParameter<Real>(
        rg.base, std::lognormal_distribution<Real>(params(0), params(1)));
  case _rdt_chi_squared:
    return RandomParameter<Real>(
        rg.base, std::chi_squared_distribution<Real>(params(0)));
  case _rdt_cauchy:
    return RandomParameter<Real>(
        rg.base, std::cauchy_distribution<Real>(params(0), params(1)));
  case _rdt_fisher_f:
    return RandomParameter<Real>(
        rg.base, std::fisher_f_distribution<Real>(params(0), params(1)));
  case _rdt_student_t:
    return RandomParameter<Real>(rg.base,
                                 std::student_t_distribution<Real>(params(0)));
  default:
    AKANTU_EXCEPTION("This is an unknown random distribution in the parser");
  }
#endif
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
