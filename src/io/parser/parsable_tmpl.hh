/**
 * @file   parsable_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Nov 13 2013
 * @date last modification: Tue Jun 24 2014
 *
 * @brief  implementation of the templated part of ParsableParam Parsable and
 * ParsableParamTyped
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

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<typename T>
const ParsableParamTyped<T> & ParsableParam::getParsableParamTyped() const {
  try {
    const ParsableParamTyped<T> & tmp = dynamic_cast<const ParsableParamTyped<T> &>(*this);
    return tmp;
  } catch (...) {
    AKANTU_EXCEPTION("The parameter named " << name << " is of type "
		     << debug::demangle(typeid(T).name()) <<".");
  }
}

/* -------------------------------------------------------------------------- */
template<typename T>
ParsableParamTyped<T> & ParsableParam::getParsableParamTyped() {
  try {
    ParsableParamTyped<T> & tmp = dynamic_cast<ParsableParamTyped<T> &>(*this);
    return tmp;
  } catch (...) {
    AKANTU_EXCEPTION("The parameter named " << name << " is of type "
		     << debug::demangle(typeid(T).name()) <<".");
  }
}

/* ------------------------------------------------------------------------ */
template<typename T, typename V>
void ParsableParam::set(const V & value) {
  ParsableParamTyped<T> & typed_param = getParsableParamTyped<T>();
  if(!(isWritable())) AKANTU_EXCEPTION("The parameter named " << name << " is not writable.");
  typed_param.setTyped(value);
}

/* -------------------------------------------------------------------------- */
template<typename T>
const T & ParsableParam::get() const {
  const ParsableParamTyped<T> & typed_param = getParsableParamTyped<T>();
  if(!(isReadable())) AKANTU_EXCEPTION("The parameter named " << name << " is not readable.");
  return typed_param.getTyped();
}

/* -------------------------------------------------------------------------- */
template<typename T>
T & ParsableParam::get() {
  ParsableParamTyped<T> & typed_param = getParsableParamTyped<T>();
  if(!(isReadable())) AKANTU_EXCEPTION("The parameter named " << name << " is not readable.");
  return typed_param.getTyped();
}


/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template<typename T>
ParsableParamTyped<T>::ParsableParamTyped(std::string name, std::string description,
					  ParamAccessType param_type, T & param) :
  ParsableParam(name, description, param_type), param(param) {}

/* -------------------------------------------------------------------------- */
template<typename T>
template<typename V>
void ParsableParamTyped<T>::setTyped(const V & value) { param = value; }

/* -------------------------------------------------------------------------- */
template<typename T>
const T & ParsableParamTyped<T>::getTyped() const { return param; }
/* -------------------------------------------------------------------------- */
template<typename T>
T & ParsableParamTyped<T>::getTyped() { return param; }

/* -------------------------------------------------------------------------- */
template<typename T>
inline void ParsableParamTyped<T>::parseParam(const ParserParameter & in_param) {
  ParsableParam::parseParam(in_param);
  param = static_cast<T>(in_param);
}

/* -------------------------------------------------------------------------- */
template<>
inline void ParsableParamTyped< Vector<Real> >::parseParam(const ParserParameter & in_param) {
  ParsableParam::parseParam(in_param);
  Vector<Real> tmp = in_param;
  for (UInt i = 0; i < param.size(); ++i) {
    param(i) = tmp(i);
  }
}

/* -------------------------------------------------------------------------- */
template<>
inline void ParsableParamTyped< Matrix<Real> >::parseParam(const ParserParameter & in_param) {
  ParsableParam::parseParam(in_param);
  Matrix<Real> tmp = in_param;
  for (UInt i = 0; i < param.cols(); ++i) {
    for (UInt j = 0; j < param.rows(); ++j) {
      param(i, j) = tmp(i, j);
    }
  }
}

/* -------------------------------------------------------------------------- */
template<>
inline void ParsableParamTyped<std::string>::parseParam(const ParserParameter & in_param) {
  ParsableParam::parseParam(in_param);
  param = in_param.getValue();
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void ParsableParamTyped<T>::printself(std::ostream & stream) const {
  ParsableParam::printself(stream);
  stream << param << std::endl;
}

/* -------------------------------------------------------------------------- */
template<>
inline void ParsableParamTyped<bool>::printself(std::ostream & stream) const {
  ParsableParam::printself(stream);
  stream << std::boolalpha << param << std::endl;
}

/* -------------------------------------------------------------------------- */
template<typename T>
void Parsable::registerParam(std::string name, T & variable,
                             ParamAccessType type,
                             const std::string description) {
  std::map<std::string, ParsableParam *>::iterator it = params.find(name);
  if(it != params.end()) AKANTU_EXCEPTION("Parameter named " << name << " already registered.");
  ParsableParamTyped<T> * param = new ParsableParamTyped<T>(name, description, type,
							    variable);
  params[name] = param;
}

/* -------------------------------------------------------------------------- */
template<typename T>
void Parsable::registerParam(std::string name, T & variable,
                             T default_value,
                             ParamAccessType type,
                             const std::string description) {
  variable = default_value;
  registerParam(name, variable, type, description);
}

/* -------------------------------------------------------------------------- */
template<typename T, typename V>
void Parsable::setMixed(const std::string & name, const V & value) {
  std::map<std::string, ParsableParam *>::iterator it = params.find(name);
  if(it == params.end()) AKANTU_EXCEPTION("No parameter named " << name << " in parsable.");
  ParsableParam & param = *(it->second);
  param.set<T>(value);
}

/* -------------------------------------------------------------------------- */
template<typename T>
void Parsable::set(const std::string & name, const T & value) {
  this->template setMixed<T>(name, value);
}

/* -------------------------------------------------------------------------- */
template<typename T>
const T & Parsable::get(const std::string & name) const {
  std::map<std::string, ParsableParam *>::const_iterator it = params.find(name);
  if(it == params.end()) AKANTU_EXCEPTION("No parameter named " << name << " in parsable.");
  const ParsableParam & param = *(it->second);
  return param.get<T>();
}

/* -------------------------------------------------------------------------- */
template<typename T>
T & Parsable::get(const std::string & name) {
  std::map<std::string, ParsableParam *>::iterator it = params.find(name);
  if(it == params.end()) AKANTU_EXCEPTION("No parameter named " << name << " in parsable.");
  ParsableParam & param = *(it->second);
  return param.get<T>();
}

__END_AKANTU__
