/**
 * @file   parser.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Nov 13 2013
 * @date last modification: Thu Aug 21 2014
 *
 * @brief  File parser interface
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
#include "aka_common.hh"
#include "aka_random_generator.hh"
/* -------------------------------------------------------------------------- */
#include <map>
/* -------------------------------------------------------------------------- */


#ifndef __AKANTU_PARSER_HH__
#define __AKANTU_PARSER_HH__

__BEGIN_AKANTU__

#define AKANTU_SECTION_TYPES                    \
  (global)                                      \
  (material)                                    \
  (model)					\
  (mesh)					\
  (heat)					\
  (contact)                                     \
  (friction)					\
  (rules)					\
  (non_local)					\
  (user)					\
  (not_defined)


#define AKANTU_SECTION_TYPES_PREFIX(elem) BOOST_PP_CAT(_st_, elem)

#define AKANTU_SECT_PREFIX(s, data, elem) AKANTU_SECTION_TYPES_PREFIX(elem)
enum SectionType {
  BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(AKANTU_SECT_PREFIX, _, AKANTU_SECTION_TYPES))
};
#undef AKANTU_SECT_PREFIX

#define AKANTU_SECTION_TYPE_PRINT_CASE(r, data, elem)                   \
  case AKANTU_SECTION_TYPES_PREFIX(elem): {                             \
    stream << BOOST_PP_STRINGIZE(AKANTU_SECTION_TYPES_PREFIX(elem));    \
    break;                                                              \
  }

inline std::ostream & operator <<(std::ostream & stream, SectionType type) {
  switch(type) {
    BOOST_PP_SEQ_FOR_EACH(AKANTU_SECTION_TYPE_PRINT_CASE, _, AKANTU_SECTION_TYPES)
  default: stream << "not a SectionType"; break;
  }
  return stream;
}

#undef AKANTU_SECTION_TYPE_PRINT_CASE
enum ParserParameterSearchCxt {
  _ppsc_current_scope             = 0x1,
  _ppsc_parent_scope              = 0x2,
  _ppsc_current_and_parent_scope  = 0x3
};

/* ------------------------------------------------------------------------ */
/* Parameters Class                                                         */
/* ------------------------------------------------------------------------ */
class ParserSection;

class ParserParameter {
public:
  ParserParameter() :
    parent_section(NULL),
    name(std::string()), value(std::string()),
    dbg_filename(std::string()) {}

  ParserParameter(const std::string & name, const std::string & value,
                  const ParserSection & parent_section) :
    parent_section(&parent_section),
    name(name), value(value),
    dbg_filename(std::string()) {}

  ParserParameter(const ParserParameter & param) :
    parent_section(param.parent_section),
    name(param.name), value(param.value),
    dbg_filename(param.dbg_filename),
    dbg_line(param.dbg_line),
    dbg_column(param.dbg_column) {}

  virtual ~ParserParameter() {}

  const std::string & getName() const { return name; }
  const std::string & getValue() const { return value; }

  void setDebugInfo(const std::string & filename,
                    UInt line, UInt column) {
    dbg_filename = filename;
    dbg_line   = line;
    dbg_column = column;
  }

  template <typename T> inline operator T() const;

  // template <typename T> inline operator Vector<T>() const;
  // template <typename T> inline operator Matrix<T>() const;


  void printself(std::ostream & stream, unsigned int indent = 0) const {
    stream << name << ": " << value
           << " (" << dbg_filename << ":" << dbg_line << ":" << dbg_column  << ")";
  }
private:
  void setParent(const ParserSection & sect) {
    parent_section = &sect;
  }

  friend class ParserSection;
private:
  const ParserSection * parent_section;
  std::string name;
  std::string value;
  std::string dbg_filename;
  UInt dbg_line, dbg_column;
};


/* ------------------------------------------------------------------------ */
/* Sections Class                                                           */
/* ------------------------------------------------------------------------ */
class ParserSection {
public:
  typedef std::multimap<SectionType, ParserSection> SubSections;
  typedef std::map<std::string, ParserParameter> Parameters;
private:
  typedef SubSections::const_iterator const_section_iterator_;
public:
  /* ------------------------------------------------------------------------ */
  /* SubSection iterator                                                      */
  /* ------------------------------------------------------------------------ */
  class const_section_iterator {
  public:
    const_section_iterator(const const_section_iterator & other) : it(other.it) { }
    const_section_iterator(const const_section_iterator_ & it) : it(it) { }

    const_section_iterator & operator=(const const_section_iterator & other) {
      if(this != &other) {
        it = other.it;
      }
      return *this;
    }
    const ParserSection & operator*  () const { return it->second; }
    const ParserSection * operator-> () const { return &(it->second); }
    bool operator==(const const_section_iterator & other) const { return it == other.it; }
    bool operator!=(const const_section_iterator & other) const { return it != other.it; }
    const_section_iterator & operator++()  { ++it; return *this; }
    const_section_iterator operator++(int) {
      const_section_iterator tmp = *this;
      operator++();
      return tmp;
    }

  private:
    const_section_iterator_ it;
  };


  /* ------------------------------------------------------------------------ */
  /* Parameters iterator                                                      */
  /* ------------------------------------------------------------------------ */
  class const_parameter_iterator {
  public:
    const_parameter_iterator(const const_parameter_iterator & other) : it(other.it) { }
    const_parameter_iterator(const Parameters::const_iterator & it) : it(it) { }

    const_parameter_iterator & operator=(const const_parameter_iterator & other) {
      if(this != &other) {
        it = other.it;
      }
      return *this;
    }
    const ParserParameter & operator* () const { return it->second; }
    const ParserParameter * operator->() { return &(it->second); };
    bool operator==(const const_parameter_iterator & other) const { return it == other.it; }
    bool operator!=(const const_parameter_iterator & other) const { return it != other.it; }
    const_parameter_iterator & operator++()  { ++it; return *this; }
    const_parameter_iterator operator++(int) { 
      const_parameter_iterator tmp = *this;
      operator++();
      return tmp;
    }

  private:
    Parameters::const_iterator it;
  };

  /* ---------------------------------------------------------------------- */
  ParserSection() :
    parent_section(NULL),
    name(std::string()), type(_st_not_defined){}

  ParserSection(const std::string & name, SectionType type) :
    parent_section(NULL), name(name), type(type) {}

  ParserSection(const std::string & name, SectionType type,
                const std::string & option,
                const ParserSection & parent_section) :
    parent_section(&parent_section), name(name), type(type), option(option) {}

  ParserSection(const ParserSection & section) :
    parent_section(section.parent_section),
    name(section.name), type(section.type),
    option(section.option),
    parameters(section.parameters),
    sub_sections_by_type(section.sub_sections_by_type) {
    setChldrenPointers();
  }

  ParserSection & operator=(const ParserSection & other) {
    if(&other != this) {
      parent_section = other.parent_section;
      name = other.name;
      type = other.type;
      option = other.option;
      parameters = other.parameters;
      sub_sections_by_type = sub_sections_by_type;
      setChldrenPointers();
    }
    return *this;
  }

  virtual ~ParserSection();

  virtual void printself(std::ostream & stream, unsigned int indent = 0) const;

  /* ---------------------------------------------------------------------- */
  /* Creation functions                                                     */
  /* ---------------------------------------------------------------------- */
public:
  ParserParameter & addParameter(const ParserParameter & param);
  ParserSection & addSubSection(const ParserSection & section);

private:
  void setChldrenPointers() {
    Parameters::iterator pit = this->parameters.begin();
    for(;pit != this->parameters.end(); ++pit) pit->second.setParent(*this);

    SubSections::iterator sit = this->sub_sections_by_type.begin();
    for (;sit != this->sub_sections_by_type.end(); ++sit) sit->second.setParent(*this);
  }

  /* ---------------------------------------------------------------------- */
  /* Accessors                                                              */
  /* ---------------------------------------------------------------------- */
public:
  std::pair<const_section_iterator, const_section_iterator>
  getSubSections(SectionType type = _st_not_defined) const {
    if(type != _st_not_defined) {
      std::pair<const_section_iterator_, const_section_iterator_> range =
	sub_sections_by_type.equal_range(type);
      return std::pair<const_section_iterator, const_section_iterator>(range.first, range.second);
    } else {
      return std::pair<const_section_iterator, const_section_iterator>(sub_sections_by_type.begin(),
								       sub_sections_by_type.end());
    }
  }

  std::pair<const_parameter_iterator, const_parameter_iterator>
  getParameters() const {
    return std::pair<const_parameter_iterator, const_parameter_iterator>(parameters.begin(),
                                                                         parameters.end());
  }

  /* ---------------------------------------------------------------------- */
  const ParserParameter & getParameter(const std::string & name,
                                       ParserParameterSearchCxt search_ctx = _ppsc_current_scope) const {
    Parameters::const_iterator it;
    if(search_ctx & _ppsc_current_scope) it = parameters.find(name);

    if(it == parameters.end()) {
      if((search_ctx & _ppsc_parent_scope) && parent_section)
        return parent_section->getParameter(name, search_ctx);
      else {
	AKANTU_SILENT_EXCEPTION("The parameter " << name
				<< " has not been found in the specified context");
      }
    }
    return it->second;
  }

  /* ------------------------------------------------------------------------ */
  bool hasParameter(const std::string & name,
		    ParserParameterSearchCxt search_ctx = _ppsc_current_scope) const {

    Parameters::const_iterator it;
    if(search_ctx & _ppsc_current_scope) it = parameters.find(name);

    if(it == parameters.end()) {
      if((search_ctx & _ppsc_parent_scope) && parent_section)
        return parent_section->hasParameter(name, search_ctx);
      else {
	return false;
      }
    }
    return true;
  }

  /* -------------------------------------------------------------------------- */
  template<class T>
  T getParameterValue(const std::string & name,
		      ParserParameterSearchCxt search_ctx = _ppsc_current_scope) const {
    const ParserParameter & tmp_param = getParameter(name, search_ctx);
    T t = tmp_param;
    return t;
  }

  /* -------------------------------------------------------------------------- */
  const std::string & getName()   const { return name; }
  const SectionType & getType()   const { return type; }
  const std::string & getOption() const { return option; }

protected:
  void setParent(const ParserSection & sect) {
    parent_section = &sect;
  }

  /* ---------------------------------------------------------------------- */
  /* Members                                                                */
  /* ---------------------------------------------------------------------- */
private:
  const ParserSection * parent_section;
  std::string name;
  SectionType type;
  std::string option;
  Parameters parameters;
  SubSections sub_sections_by_type;
};


/* ------------------------------------------------------------------------ */
/* Parser Class                                                             */
/* ------------------------------------------------------------------------ */
class Parser : public ParserSection {
public:
  Parser() : ParserSection("global", _st_global) {}

  void parse(const std::string & filename);

  std::string getLastParsedFile() const;

  static bool isPermissive() { return parser_permissive; }
public:

  static Real         parseReal  (const std::string & value, const ParserSection & section);
  static Vector<Real> parseVector(const std::string & value, const ParserSection & section);
  static Matrix<Real> parseMatrix(const std::string & value, const ParserSection & section);
  static RandomParameter<Real> parseRandomParameter(const std::string & value, const ParserSection & section);
protected:
  template <class T, class Grammar>
  static T parseType(const std::string & value, Grammar & grammar);

protected:
  friend class Parsable;
  static bool parser_permissive;
  std::string last_parsed_file;
};

inline std::ostream & operator <<(std::ostream & stream, const ParserParameter &_this) {
  _this.printself(stream);
  return stream;
}

inline std::ostream & operator <<(std::ostream & stream, const ParserSection & section) {
  section.printself(stream);
  return stream;
}

__END_AKANTU__


#include "parser_tmpl.hh"


#endif /* __AKANTU_PARSER_HH__ */
