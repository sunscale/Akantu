/**
 * @file   parameter_registry.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed May 04 2016
 * @date last modification: Thu Feb 01 2018
 *
 * @brief  Parameter Registry and derived classes implementation
 *
 * @section LICENSE
 *
 * Copyright (©) 2016-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include <utility>

#include "parameter_registry.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

Parameter::Parameter() : name(""), description("") {}

/* -------------------------------------------------------------------------- */
Parameter::Parameter(std::string name, std::string description,
                     ParameterAccessType param_type)
    : name(std::move(name)), description(std::move(description)),
      param_type(param_type) {}

/* -------------------------------------------------------------------------- */
bool Parameter::isWritable() const { return param_type & _pat_writable; }

/* -------------------------------------------------------------------------- */
bool Parameter::isReadable() const { return param_type & _pat_readable; }

/* -------------------------------------------------------------------------- */
bool Parameter::isInternal() const { return param_type & _pat_internal; }

/* -------------------------------------------------------------------------- */
bool Parameter::isParsable() const { return param_type & _pat_parsable; }

/* -------------------------------------------------------------------------- */
void Parameter::setAccessType(ParameterAccessType ptype) {
  this->param_type = ptype;
}

/* -------------------------------------------------------------------------- */
void Parameter::printself(std::ostream & stream) const {
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
/* -------------------------------------------------------------------------- */
ParameterRegistry::ParameterRegistry() = default;

/* -------------------------------------------------------------------------- */
ParameterRegistry::~ParameterRegistry() {
  std::map<std::string, Parameter *>::iterator it, end;
  for (it = params.begin(); it != params.end(); ++it) {
    delete it->second;
    it->second = NULL;
  }
  this->params.clear();
}

/* -------------------------------------------------------------------------- */
void ParameterRegistry::printself(std::ostream & stream, int indent) const {
  std::string space(indent, AKANTU_INDENT);

  Parameters::const_iterator it;
  for (it = params.begin(); it != params.end(); ++it) {
    stream << space;
    it->second->printself(stream);
  }

  SubRegisteries::const_iterator sub_it;
  for (sub_it = sub_registries.begin(); sub_it != sub_registries.end();
       ++sub_it) {
    stream << space << "Registry [" << std::endl;
    sub_it->second->printself(stream, indent + 1);
    stream << space << "]";
  }
}

/* -------------------------------------------------------------------------- */
void ParameterRegistry::registerSubRegistry(const ID & id,
                                            ParameterRegistry & registry) {
  sub_registries[id] = &registry;
}

/* -------------------------------------------------------------------------- */
void ParameterRegistry::setParameterAccessType(const std::string & name,
                                               ParameterAccessType ptype) {
  auto it = params.find(name);
  if (it == params.end())
    AKANTU_CUSTOM_EXCEPTION(debug::ParameterUnexistingException(name, *this));
  Parameter & param = *(it->second);
  param.setAccessType(ptype);
}

} // namespace akantu
