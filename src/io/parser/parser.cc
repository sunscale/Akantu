/**
 * @file   parser.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Nov 13 2013
 * @date last modification: Thu Feb 01 2018
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
// STL
#include <fstream>
#include <iomanip>
#include <map>

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "parser.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
ParserSection::~ParserSection() { this->clean(); }

/* -------------------------------------------------------------------------- */
ParserParameter & ParserSection::addParameter(const ParserParameter & param) {
  if (parameters.find(param.getName()) != parameters.end())
    AKANTU_EXCEPTION("The parameter \"" + param.getName() +
                     "\" is already defined in this section");

  return (parameters
              .insert(std::pair<std::string, ParserParameter>(param.getName(),
                                                              param))
              .first->second);
}

/* -------------------------------------------------------------------------- */
ParserSection & ParserSection::addSubSection(const ParserSection & section) {
  return ((sub_sections_by_type.insert(std::pair<ParserType, ParserSection>(
               section.getType(), section)))
              ->second);
}

/* -------------------------------------------------------------------------- */
std::string Parser::getLastParsedFile() const { return last_parsed_file; }

/* -------------------------------------------------------------------------- */
void ParserSection::printself(std::ostream & stream,
                              unsigned int indent) const {
  std::string space(indent, AKANTU_INDENT);
  stream << space << "Section(" << this->type << ") " << this->name
         << (option != "" ? (" " + option) : "") << " [" << std::endl;
  if (!this->parameters.empty()) {
    stream << space << " Parameters [" << std::endl;
    auto pit = this->parameters.begin();
    for (; pit != this->parameters.end(); ++pit) {
      stream << space << " + ";
      pit->second.printself(stream, indent);
      stream << "\n";
    }
    stream << space << " ]" << std::endl;
  }

  if (!this->sub_sections_by_type.empty()) {
    stream << space << " Subsections [" << std::endl;
    auto sit = this->sub_sections_by_type.begin();
    for (; sit != this->sub_sections_by_type.end(); ++sit)
      sit->second.printself(stream, indent + 2);
    stream << std::endl;
    stream << space << " ]" << std::endl;
  }
  stream << space << "]" << std::endl;
}

} // namespace akantu
