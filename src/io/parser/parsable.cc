/**
 * @file   parsable.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Nov 13 2013
 * @date last modification: Thu Nov 19 2015
 *
 * @brief  Parsable implementation
 *
 * @section LICENSE
 *
 * Copyright  (©)  2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de Lausanne)
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
#include "parsable.hh"
#include "aka_random_generator.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

ParsableParam::ParsableParam()
    : name(""), description(""), param_type(_pat_internal) {}

/* -------------------------------------------------------------------------- */
ParsableParam::ParsableParam(std::string name, std::string description,
                             ParamAccessType param_type)
    : name(name), description(description), param_type(param_type) {}

/* -------------------------------------------------------------------------- */
bool ParsableParam::isWritable() const { return param_type & _pat_writable; }

/* -------------------------------------------------------------------------- */
bool ParsableParam::isReadable() const { return param_type & _pat_readable; }

/* -------------------------------------------------------------------------- */
bool ParsableParam::isInternal() const { return param_type & _pat_internal; }

/* -------------------------------------------------------------------------- */
bool ParsableParam::isParsable() const { return param_type & _pat_parsable; }

/* -------------------------------------------------------------------------- */
void ParsableParam::setAccessType(ParamAccessType ptype) {
  this->param_type = ptype;
}

/* -------------------------------------------------------------------------- */
void ParsableParam::printself(std::ostream & stream) const {
  stream << " ";
  if (isInternal())
    stream << "iii";
  else {
    if (isReadable())
      stream << "r";
    else
      stream << "-";

    if (isWritable())
      stream << "w";
    else
      stream << "-";

    if (isParsable())
      stream << "p";
    else
      stream << "-";
  }
  stream << " ";

  std::stringstream sstr;
  sstr << name;
  UInt width = std::max(int(10 - sstr.str().length()), 0);
  sstr.width(width);

  if (description != "") {
    sstr << " [" << description << "]";
  }

  stream << sstr.str();
  width = std::max(int(50 - sstr.str().length()), 0);
  stream.width(width);

  stream << " : ";
}

/* -------------------------------------------------------------------------- */
void ParsableParam::parseParam(const ParserParameter & param) {
  if (!isParsable())
    AKANTU_EXCEPTION("The parameter named " << param.getName()
                                            << " is not parsable.");
}

/* -------------------------------------------------------------------------- */
Parsable::~Parsable() {
  std::map<std::string, ParsableParam *>::iterator it, end;
  for (it = params.begin(); it != params.end(); ++it) {
    delete it->second;
    it->second = NULL;
  }
  this->params.clear();
}

/* -------------------------------------------------------------------------- */
void Parsable::printself(std::ostream & stream, int indent) const {
  std::string space;
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
    ;
  std::map<std::string, ParsableParam *>::const_iterator it, end;
  for (it = params.begin(); it != params.end(); ++it) {
    stream << space;
    it->second->printself(stream);
  }

  SubSections::const_iterator sit, send;
  for (sit = sub_sections.begin(); sit != sub_sections.end(); ++sit) {
    sit->second->printself(stream, indent + 1);
  }
}

/* -------------------------------------------------------------------------- */
void Parsable::setParamAccessType(const std::string & name,
                                  ParamAccessType ptype) {
  std::map<std::string, ParsableParam *>::iterator it = params.find(name);
  if (it == params.end())
    AKANTU_EXCEPTION("No parameter named " << name << " in parsable.");
  ParsableParam & param = *(it->second);
  param.setAccessType(ptype);
}

/* -------------------------------------------------------------------------- */
void Parsable::registerSubSection(const SectionType & type,
                                  const std::string & name,
                                  Parsable & sub_section) {
  SubSectionKey key(type, name);
  sub_sections[key] = &sub_section;
}

/* -------------------------------------------------------------------------- */
void Parsable::parseParam(const ParserParameter & in_param) {
  std::map<std::string, ParsableParam *>::iterator it =
      params.find(in_param.getName());
  if (it == params.end()) {
    if (Parser::isPermissive()) {
      AKANTU_DEBUG_WARNING("No parameter named " << in_param.getName()
                                                 << " registered in " << pid
                                                 << ".");
      return;
    } else
      AKANTU_EXCEPTION("No parameter named " << in_param.getName()
                                             << " registered in " << pid
                                             << ".");
  }
  ParsableParam & param = *(it->second);
  param.parseParam(in_param);
}

/* -------------------------------------------------------------------------- */
void Parsable::parseSection(const ParserSection & section) {
  if (section_type != section.getType())
    AKANTU_EXCEPTION("The object "
                     << pid << " is meant to parse section of type "
                     << section_type << ", so it cannot parse section of type "
                     << section.getType());

  std::pair<Parser::const_parameter_iterator, Parser::const_parameter_iterator>
      params = section.getParameters();
  Parser::const_parameter_iterator it = params.first;
  for (; it != params.second; ++it) {
    parseParam(*it);
  }

  Parser::const_section_iterator sit = section.getSubSections().first;
  for (; sit != section.getSubSections().second; ++sit) {
    parseSubSection(*sit);
  }
}

/* -------------------------------------------------------------------------- */
void Parsable::parseSubSection(const ParserSection & section) {
  SubSectionKey key(section.getType(), section.getName());
  SubSections::iterator it = sub_sections.find(key);
  if (it != sub_sections.end()) {
    it->second->parseSection(section);
  } else if (!Parser::isPermissive()) {
    AKANTU_EXCEPTION("No parsable defined for sub sections of type <"
                     << key.first << "," << key.second << "> in " << pid);
  } else
    AKANTU_DEBUG_WARNING("No parsable defined for sub sections of type <"
                         << key.first << "," << key.second << "> in " << pid);
}

__END_AKANTU__
