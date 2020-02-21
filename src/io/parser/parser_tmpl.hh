/**
 * @file   parser_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Nov 13 2013
 * @date last modification: Mon Dec 18 2017
 *
 * @brief  Implementation of the parser templated methods
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
#include <regex>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <typename T> inline ParserParameter::operator T() const {
  T t;
  std::stringstream sstr(value);
  sstr >> t;
  if (sstr.bad())
    AKANTU_EXCEPTION("No known conversion of a ParserParameter \""
                     << name << "\" to the type " << typeid(T).name());
  return t;
}

/* -------------------------------------------------------------------------- */
template <> inline ParserParameter::operator const char *() const {
  return value.c_str();
}

/* -------------------------------------------------------------------------- */
template <> inline ParserParameter::operator Real() const {
  return Parser::parseReal(value, *parent_section);
}

/* --------------------------------------------------------- -----------------
 */
template <> inline ParserParameter::operator bool() const {
  bool b;
  std::stringstream sstr(value);
  sstr >> std::boolalpha >> b;
  if (sstr.fail()) {
    sstr.clear();
    sstr >> std::noboolalpha >> b;
  }
  return b;
}

/* -------------------------------------------------------------------------- */
template <> inline ParserParameter::operator std::vector<std::string>() const {
  std::vector<std::string> tmp;
  auto string =
      std::regex_replace(value, std::regex("[[:space:]]|\\[|\\]"), "");
  std::smatch sm;
  while (std::regex_search(string, sm, std::regex("[^,]+"))) {
    tmp.push_back(sm.str());
    string = sm.suffix();
  }
  return tmp;
}

/* -------------------------------------------------------------------------- */
template <> inline ParserParameter::operator std::set<std::string>() const {
  std::set<std::string> tmp;
  auto string =
      std::regex_replace(value, std::regex("[[:space:]]|\\[|\\]"), "");
  std::smatch sm;
  while (std::regex_search(string, sm, std::regex("[^,]+"))) {
    tmp.emplace(sm.str());
    string = sm.suffix();
  }
  return tmp;
}

/* -------------------------------------------------------------------------- */
template <> inline ParserParameter::operator Vector<Real>() const {
  return Parser::parseVector(value, *parent_section);
}

/* --------------------------------------------------------- -----------------
 */
template <> inline ParserParameter::operator Vector<UInt>() const {
  Vector<Real> tmp = Parser::parseVector(value, *parent_section);
  Vector<UInt> tmp_uint(tmp.size());
  for (UInt i = 0; i < tmp.size(); ++i) {
    tmp_uint(i) = UInt(tmp(i));
  }
  return tmp_uint;
}

/* --------------------------------------------------------- -----------------
 */
template <> inline ParserParameter::operator Matrix<Real>() const {
  return Parser::parseMatrix(value, *parent_section);
}

/* -------------------------------------------------------------------------- */
template <> inline ParserParameter::operator RandomParameter<Real>() const {
  return Parser::parseRandomParameter(value, *parent_section);
}

} // namespace akantu
