/**
 * @file   parsable.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Aug 09 2012
 * @date last modification: Fri Dec 08 2017
 *
 * @brief  Interface of the parameter registry
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_common.hh"
#include "parameter_registry.hh"
#include "parser.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_PARSABLE_HH_
#define AKANTU_PARSABLE_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
/* Parsable Interface                                                         */
/* -------------------------------------------------------------------------- */
/// Defines interface for classes to manipulate parsable parameters
class Parsable : public ParameterRegistry {
public:
  Parsable(const ParserType & section_type, const ID & id = std::string());
  ~Parsable() override;

  /// Add subsection to the sub_sections map
  void registerSubSection(const ParserType & type, const std::string & name,
                          Parsable & sub_section);

  /* ------------------------------------------------------------------------ */
public:
  virtual void parseSection(const ParserSection & section);
  virtual void parseSubSection(const ParserSection & section);
  virtual void parseParam(const ParserParameter & in_param);

private:
  ParserType section_type;
  /// ID of parsable object
  ID pid;
  using SubSectionKey = std::pair<ParserType, std::string>;
  using SubSections = std::map<SubSectionKey, Parsable *>;
  /// Subsections map
  SubSections sub_sections;
};

} // namespace akantu

#endif /* AKANTU_PARSABLE_HH_ */
