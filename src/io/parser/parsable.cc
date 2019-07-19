/**
 * @file   parsable.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Nov 13 2013
 * @date last modification: Thu Feb 08 2018
 *
 * @brief  Parsable implementation
 *
 * @section LICENSE
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
#include "parsable.hh"
#include "aka_random_generator.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
Parsable::Parsable(const ParserType & section_type, const ID & id)
    : section_type(section_type), pid(id) {
  this->consisder_sub = false;
}

/* -------------------------------------------------------------------------- */
Parsable::~Parsable() = default;

/* -------------------------------------------------------------------------- */
void Parsable::registerSubSection(const ParserType & type,
                                  const std::string & name,
                                  Parsable & sub_section) {
  SubSectionKey key(type, name);
  sub_sections[key] = &sub_section;

  this->registerSubRegistry(name, sub_section);
}

/* -------------------------------------------------------------------------- */
void Parsable::parseParam(const ParserParameter & in_param) {
  auto it = params.find(in_param.getName());
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
  Parameter & param = *(it->second);
  param.setAuto(in_param);
}

/* -------------------------------------------------------------------------- */
void Parsable::parseSection(const ParserSection & section) {
  if (section_type != section.getType())
    AKANTU_EXCEPTION("The object "
                     << pid << " is meant to parse section of type "
                     << section_type << ", so it cannot parse section of type "
                     << section.getType());

  auto params = section.getParameters();
  auto it = params.first;
  for (; it != params.second; ++it) {
    parseParam(*it);
  }

  auto sit = section.getSubSections().first;
  for (; sit != section.getSubSections().second; ++sit) {
    parseSubSection(*sit);
  }
}

/* -------------------------------------------------------------------------- */
void Parsable::parseSubSection(const ParserSection & section) {
  SubSectionKey key(section.getType(), section.getName());
  auto it = sub_sections.find(key);
  if (it != sub_sections.end()) {
    it->second->parseSection(section);
  } else if (!Parser::isPermissive()) {
    AKANTU_EXCEPTION("No parsable defined for sub sections of type <"
                     << key.first << "," << key.second << "> in " << pid);
  } else {
    AKANTU_DEBUG_WARNING("No parsable defined for sub sections of type <"
                         << key.first << "," << key.second << "> in " << pid);
  }
}

} // namespace akantu
