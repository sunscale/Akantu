/**
 * @file   parsable.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Aug 09 2012
 * @date last modification: Thu Dec 17 2015
 *
 * @brief  Interface of parsable objects
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "aka_common.hh"
#include "parser.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_PARSABLE_HH__
#define __AKANTU_PARSABLE_HH__

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/// Defines the access modes of parsable parameters
enum ParamAccessType {
  _pat_internal   = 0x0001,
  _pat_writable   = 0x0010,
  _pat_readable   = 0x0100,
  _pat_modifiable = 0x0110, //_pat_readable | _pat_writable,
  _pat_parsable   = 0x1000,
  _pat_parsmod    = 0x1110  //< _pat_parsable | _pat_modifiable
};


/// Bit-wise operator between access modes
inline ParamAccessType operator|(const ParamAccessType & a, const ParamAccessType & b) {
  ParamAccessType tmp = ParamAccessType(UInt(a) | UInt(b));
  return tmp;
}
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template<typename T> class ParsableParamTyped;

/* -------------------------------------------------------------------------- */
/**
 * Interface for the ParsableParamTyped
 */
class ParsableParam {
public:
  ParsableParam();
  ParsableParam(std::string name, std::string description, ParamAccessType param_type);

  virtual ~ParsableParam() {};
  /* ------------------------------------------------------------------------ */
  bool isInternal() const;
  bool isWritable() const;
  bool isReadable() const;
  bool isParsable() const;

  void setAccessType(ParamAccessType ptype);

  /* ------------------------------------------------------------------------ */
  template<typename T, typename V> void set(const V & value);
  template<typename T> T & get();
  template<typename T> const T & get() const;

  virtual void parseParam(const ParserParameter & param);

  /* ------------------------------------------------------------------------ */
  virtual void printself(std::ostream & stream) const;

protected:
  /// Returns const instance of templated sub-class ParsableParamTyped
  template<typename T>
  const ParsableParamTyped<T> & getParsableParamTyped() const;

  /// Returns instance of templated sub-class ParsableParamTyped
  template<typename T>
  ParsableParamTyped<T> & getParsableParamTyped();

private:
  /// Name of parameter
  std::string name;
  /// Description of parameter
  std::string description;
  /// Type of access
  ParamAccessType param_type;
};

/* -------------------------------------------------------------------------- */
/* Typed Parameter                                                            */
/* -------------------------------------------------------------------------- */
/**
 * Type parameter transfering a ParserParameter (string: string) to a typed parameter in the memory of the p
 */
template<typename T>
class ParsableParamTyped : public ParsableParam {
public:
  ParsableParamTyped(std::string name, std::string description,
		     ParamAccessType param_type, T & param);

  /* ------------------------------------------------------------------------ */
  template<typename V>
  void setTyped(const V & value);
  T & getTyped();
  const T & getTyped() const;

  void parseParam(const ParserParameter & param);

  virtual void printself(std::ostream & stream) const;
private:
  /// Value of parameter
  T & param;
};


/* -------------------------------------------------------------------------- */
/* Parsable Interface                                                         */
/* -------------------------------------------------------------------------- */
/// Defines interface for classes to manipulate parsable parameters
class Parsable {
public:
  Parsable(const SectionType & section_type,
           const ID & id = std::string()) : section_type(section_type), pid(id) {};
  virtual ~Parsable();

  /* ------------------------------------------------------------------------ */
  /// Add parameter to the params map
  template<typename T>
  void registerParam(std::string name, T & variable,
		     ParamAccessType type,
		     const std::string description = "");

  /// Add parameter to the params map (with default value)
  template<typename T>
  void registerParam(std::string name, T & variable, T default_value,
		     ParamAccessType type,
		     const std::string description = "");

  /// Add subsection to the sub_sections map
  void registerSubSection(const SectionType & type,
			  const std::string & name,
			  Parsable & sub_section);

  /* ------------------------------------------------------------------------ */
  /// Set value to a parameter (with possible different type)
  template<typename T, typename V> void setMixed(const std::string & name, const V & value);
  /// Set value to a parameter
  template<typename T> void set(const std::string & name, const T & value);
  /// Get value of a parameter
  template<typename T> const T & get(const std::string & name) const;
protected:
  template<typename T> T & get(const std::string & name);

protected:
  void setParamAccessType(const std::string & name, ParamAccessType ptype);
  /* ------------------------------------------------------------------------ */
public:
  virtual void parseSection(const ParserSection & section);

  virtual void parseSubSection(const ParserSection & section);

  virtual void parseParam(const ParserParameter & parameter);
  /* ------------------------------------------------------------------------ */
  virtual void printself(std::ostream & stream, int indent) const;

private:
  SectionType section_type;
  /// ID of parsable object
  ID pid;
  /// Parameters map
  std::map<std::string, ParsableParam *> params;
  typedef std::pair<SectionType, std::string> SubSectionKey;
  typedef std::map<SubSectionKey, Parsable *> SubSections;
  /// Subsections map
  SubSections sub_sections;
};

__END_AKANTU__

#include "parsable_tmpl.hh"


#endif /* __AKANTU_PARSABLE_HH__ */
