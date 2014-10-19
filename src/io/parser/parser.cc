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
// STL
#include <fstream>
#include <iomanip>
#include <map>

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "parser.hh"
#include "static_communicator.hh"
/* -------------------------------------------------------------------------- */
//#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/support_multi_pass.hpp>
#include <boost/spirit/include/classic_position_iterator.hpp>
/* -------------------------------------------------------------------------- */
#include "algebraic_parser.hh"
#include "input_file_parser.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
ParserSection::~ParserSection() {
}

/* -------------------------------------------------------------------------- */
ParserParameter & ParserSection::addParameter(const ParserParameter & param) {
  if(parameters.find(param.getName()) != parameters.end())
    AKANTU_EXCEPTION("The parameter \"" + param.getName() + "\" is already defined in this section");

  return (parameters.insert
	  (std::pair<std::string, ParserParameter>(param.getName(),
						   param)).first->second);
}

/* -------------------------------------------------------------------------- */
ParserSection & ParserSection::addSubSection(const ParserSection & section) {
  return ((sub_sections_by_type.insert
	   (std::pair<SectionType, ParserSection>(section.getType(),
						  section)))->second);
}

/* -------------------------------------------------------------------------- */
std::string Parser::getLastParsedFile() const {
  return last_parsed_file;
}


/* -------------------------------------------------------------------------- */
void Parser::parse(const std::string& filename) {
  std::ifstream input(filename.c_str());

  if(!input.good()) {
    AKANTU_EXCEPTION("Could not open file " << filename << "!");
  }

  input.unsetf(std::ios::skipws);

  // wrap istream into iterator
  spirit::istream_iterator fwd_begin(input);
  spirit::istream_iterator fwd_end;

  // wrap forward iterator with position iterator, to record the position
  typedef spirit::classic::position_iterator2<spirit::istream_iterator> pos_iterator_type;
  pos_iterator_type position_begin(fwd_begin, fwd_end, filename);
  pos_iterator_type position_end;

  // parse
  parser::InputFileGrammar<pos_iterator_type> ag(this);

  bool result = qi::phrase_parse(position_begin,
				 position_end,
				 ag,
				 ag.skipper);

  if(!result || position_begin != position_end) {
    spirit::classic::file_position pos = position_begin.get_position();

    AKANTU_EXCEPTION("Parse error [ " << ag.getErrorMessage() << " ]"
		     << " in file " << filename
		     << " line " << pos.line
		     << " column " << pos.column << std::endl
		     << "'" << position_begin.get_currentline() << "'" << std::endl
		     << std::setw(pos.column) << " " << "^- here");
  }


  try {
    DebugLevel dbl = debug::getDebugLevel();
    debug::setDebugLevel(dblError);
    bool permissive = getParameter("parser_permissive", _ppsc_current_scope);
    debug::setDebugLevel(dbl);

    parser_permissive = permissive;
    AKANTU_DEBUG_INFO("Parser switched permissive mode to " << std::boolalpha << parser_permissive);
  } catch (debug::Exception & e) {
  }

  last_parsed_file = filename;
  input.close();
}

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


/* -------------------------------------------------------------------------- */
Real Parser::parseReal(const std::string & value, const ParserSection & section) {
  using boost::spirit::ascii::space_type;
  parser::AlgebraicGrammar<std::string::const_iterator, space_type> grammar(section);
  return Parser::parseType<Real>(value, grammar);
}

/* -------------------------------------------------------------------------- */
Vector<Real> Parser::parseVector(const std::string & value, const ParserSection & section) {
  using boost::spirit::ascii::space_type;
  parser::VectorGrammar<std::string::const_iterator, space_type> grammar(section);
  return Parser::parseType<parser::parsable_vector>(value, grammar);
}

/* -------------------------------------------------------------------------- */
Matrix<Real> Parser::parseMatrix(const std::string & value, const ParserSection & section) {
  using boost::spirit::ascii::space_type;
  parser::MatrixGrammar<std::string::const_iterator, space_type> grammar(section);
  return Parser::parseType<parser::parsable_matrix>(value, grammar);
}

/* -------------------------------------------------------------------------- */
RandomParameter<Real> Parser::parseRandomParameter(const std::string & value, const ParserSection & section) {
  using boost::spirit::ascii::space_type;
  parser::RandomGeneratorGrammar<std::string::const_iterator, space_type> grammar(section);
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
void ParserSection::printself(std::ostream & stream, unsigned int indent) const {
  std::string space;
  std::string ind = AKANTU_INDENT;
  for(unsigned int i = 0; i < indent; i++, space += ind);

  stream << space << "Section(" << this->type << ") "
	 << this->name << (option != "" ? (" " + option) : "")
	 << " [" << std::endl;
  if(!this->parameters.empty()) {
    stream << space << ind << "Parameters [" << std::endl;
    Parameters::const_iterator pit = this->parameters.begin();
    for(;pit != this->parameters.end(); ++pit) {
      stream << space << ind << " + ";
      pit->second.printself(stream);
      stream << std::endl;
    }
    stream << space << ind << "]" << std::endl;
  }

  if(!this->sub_sections_by_type.empty()) {
    stream << space << ind << "Subsections [" << std::endl;
    SubSections::const_iterator sit = this->sub_sections_by_type.begin();
    for (;sit != this->sub_sections_by_type.end(); ++sit) sit->second.printself(stream, indent + 2);
    stream << std::endl;
    stream << space << ind << "]" << std::endl;
  }
  stream << space << "]" << std::endl;
}

__END_AKANTU__
