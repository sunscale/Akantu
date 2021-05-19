/**
 * @file   parser.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Nov 13 2013
 * @date last modification: Fri Dec 08 2017
 *
 * @brief  File parser interface
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
#include "aka_common.hh"
#include "aka_random_generator.hh"
/* -------------------------------------------------------------------------- */
#include <map>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_PARSER_HH_
#define AKANTU_PARSER_HH_

namespace akantu {

#if !defined(DOXYGEN)
// clang-format off
#define AKANTU_SECTION_TYPES                                            \
  (cohesive_inserter)                                                   \
  (contact)                                                             \
  (embedded_interface)                                                  \
  (friction)                                                            \
  (global)                                                              \
  (heat)                                                                \
  (integration_scheme)                                                  \
  (material)                                                            \
  (phasefield)								\
  (mesh)                                                                \
  (model)                                                               \
  (model_solver)                                                        \
  (neighborhood)                                                        \
  (neighborhoods)                                                       \
  (non_linear_solver)                                                   \
  (non_local)                                                           \
  (rules)                                                               \
  (solver)                                                              \
  (time_step_solver)                                                    \
  (user)                                                                \
  (weight_function)                                                     \
  (not_defined)
// clang-format on

/// Defines the possible section types
AKANTU_CLASS_ENUM_DECLARE(ParserType, AKANTU_SECTION_TYPES)
AKANTU_CLASS_ENUM_OUTPUT_STREAM(ParserType, AKANTU_SECTION_TYPES)
AKANTU_CLASS_ENUM_INPUT_STREAM(ParserType, AKANTU_SECTION_TYPES)
#endif

/// Defines the possible search contexts/scopes (for parameter search)
enum ParserParameterSearchCxt {
  _ppsc_current_scope = 0x1,
  _ppsc_parent_scope = 0x2,
  _ppsc_current_and_parent_scope = 0x3
};
/* ------------------------------------------------------------------------ */
/* Parameters Class                                                         */
/* ------------------------------------------------------------------------ */
class ParserSection;

/// @brief The ParserParameter objects represent the end of tree branches as
/// they
/// are the different informations contained in the input file.
class ParserParameter {
public:
  ParserParameter()
      : name(std::string()), value(std::string()), dbg_filename(std::string()) {
  }

  ParserParameter(const std::string & name, const std::string & value,
                  const ParserSection & parent_section)
      : parent_section(&parent_section), name(name), value(value),
        dbg_filename(std::string()) {}

  ParserParameter(const ParserParameter & param) = default;

  virtual ~ParserParameter() = default;

  /// Get parameter name
  const std::string & getName() const { return name; }
  /// Get parameter value
  const std::string & getValue() const { return value; }

  /// Set info for debug output
  void setDebugInfo(const std::string & filename, UInt line, UInt column) {
    dbg_filename = filename;
    dbg_line = line;
    dbg_column = column;
  }

  template <typename T> inline operator T() const;

  // template <typename T> inline operator Vector<T>() const;
  // template <typename T> inline operator Matrix<T>() const;

  /// Print parameter info in stream
  void printself(std::ostream & stream,
                 __attribute__((unused)) unsigned int indent = 0) const {
    stream << name << ": " << value << " (" << dbg_filename << ":" << dbg_line
           << ":" << dbg_column << ")";
  }

private:
  void setParent(const ParserSection & sect) { parent_section = &sect; }

  friend class ParserSection;

private:
  /// Pointer to the parent section
  const ParserSection * parent_section{nullptr};
  /// Name of the parameter
  std::string name;
  /// Value of the parameter
  std::string value;
  /// File for debug output
  std::string dbg_filename;
  /// Position of parameter in parsed file
  UInt dbg_line, dbg_column;
};

/* ------------------------------------------------------------------------ */
/* Sections Class                                                           */
/* ------------------------------------------------------------------------ */
/// ParserSection represents a branch of the parsing tree.
class ParserSection {
public:
  using SubSections = std::multimap<ParserType, ParserSection>;
  using Parameters = std::map<std::string, ParserParameter>;

private:
  using const_section_iterator_ = SubSections::const_iterator;

public:
  /* ------------------------------------------------------------------------ */
  /* SubSection iterator                                                      */
  /* ------------------------------------------------------------------------ */
  /// Iterator on sections
  class const_section_iterator {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = ParserSection;
    using pointer = ParserSection *;
    using reference = ParserSection &;

    const_section_iterator() = default;
    const_section_iterator(const const_section_iterator_ & it) : it(it) {}
    const_section_iterator(const const_section_iterator & other) = default;
    const_section_iterator &
    operator=(const const_section_iterator & other) = default;

    const ParserSection & operator*() const { return it->second; }
    const ParserSection * operator->() const { return &(it->second); }

    bool operator==(const const_section_iterator & other) const {
      return it == other.it;
    }

    bool operator!=(const const_section_iterator & other) const {
      return it != other.it;
    }

    const_section_iterator & operator++() {
      ++it;
      return *this;
    }

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
  /// Iterator on parameters
  class const_parameter_iterator {
  public:
    const_parameter_iterator(const const_parameter_iterator & other) = default;
    const_parameter_iterator(const Parameters::const_iterator & it) : it(it) {}

    const_parameter_iterator &
    operator=(const const_parameter_iterator & other) {
      if (this != &other) {
        it = other.it;
      }
      return *this;
    }

    const ParserParameter & operator*() const { return it->second; }
    const ParserParameter * operator->() { return &(it->second); };

    bool operator==(const const_parameter_iterator & other) const {
      return it == other.it;
    }

    bool operator!=(const const_parameter_iterator & other) const {
      return it != other.it;
    }

    const_parameter_iterator & operator++() {
      ++it;
      return *this;
    }

    const_parameter_iterator operator++(int) {
      const_parameter_iterator tmp = *this;
      operator++();
      return tmp;
    }

  private:
    Parameters::const_iterator it;
  };

  /* ---------------------------------------------------------------------- */
  ParserSection() : name(std::string()) {}

  ParserSection(const std::string & name, ParserType type)
      : name(name), type(type) {}

  ParserSection(const std::string & name, ParserType type,
                const std::string & option,
                const ParserSection & parent_section)
      : parent_section(&parent_section), name(name), type(type),
        option(option) {}

  ParserSection(const ParserSection & section)
      : parent_section(section.parent_section), name(section.name),
        type(section.type), option(section.option),
        parameters(section.parameters),
        sub_sections_by_type(section.sub_sections_by_type) {
    setChldrenPointers();
  }

  ParserSection & operator=(const ParserSection & other) {
    if (&other != this) {
      parent_section = other.parent_section;
      name = other.name;
      type = other.type;
      option = other.option;
      parameters = other.parameters;
      sub_sections_by_type = other.sub_sections_by_type;
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

protected:
  /// Clean ParserSection content
  void clean() {
    parameters.clear();
    sub_sections_by_type.clear();
  }

private:
  void setChldrenPointers() {
    for (auto && param_pair : this->parameters) {
      param_pair.second.setParent(*this);
    }

    for (auto && sub_sect_pair : this->sub_sections_by_type) {
      sub_sect_pair.second.setParent(*this);
    }
  }

  /* ---------------------------------------------------------------------- */
  /* Accessors                                                              */
  /* ---------------------------------------------------------------------- */
public:
  class SubSectionsRange
      : public std::pair<const_section_iterator, const_section_iterator> {
  public:
    SubSectionsRange(const const_section_iterator & first,
                     const const_section_iterator & second)
        : std::pair<const_section_iterator, const_section_iterator>(first,
                                                                    second) {}
    auto begin() { return this->first; }
    auto end() { return this->second; }
  };

  /// Get begin and end iterators on subsections of certain type
  auto getSubSections(ParserType type = ParserType::_not_defined) const {
    if (type != ParserType::_not_defined) {
      auto range = sub_sections_by_type.equal_range(type);
      return SubSectionsRange(range.first, range.second);
    }
    return SubSectionsRange(sub_sections_by_type.begin(),
                            sub_sections_by_type.end());
  }

  /// Get number of subsections of certain type
  UInt getNbSubSections(ParserType type = ParserType::_not_defined) const {
    if (type != ParserType::_not_defined) {
      return this->sub_sections_by_type.count(type);
    }
    return this->sub_sections_by_type.size();
  }

  /// Get begin and end iterators on parameters
  auto getParameters() const {
    return std::pair<const_parameter_iterator, const_parameter_iterator>(
        parameters.begin(), parameters.end());
  }
  /* ---------------------------------------------------------------------- */
  /// Get parameter within specified context
  const ParserParameter & getParameter(
      const std::string & name,
      ParserParameterSearchCxt search_ctx = _ppsc_current_scope) const {
    Parameters::const_iterator it;
    if ((search_ctx & _ppsc_current_scope) != 0) {
      it = parameters.find(name);
    }

    if (it == parameters.end()) {
      if ((search_ctx & _ppsc_parent_scope) != 0 and
          parent_section != nullptr) {
        return parent_section->getParameter(name, search_ctx);
      }

      AKANTU_SILENT_EXCEPTION(
          "The parameter " << name
                           << " has not been found in the specified context");
    }
    return it->second;
  }

  /* ------------------------------------------------------------------------ */
  /// Get parameter within specified context, with a default value in case the
  /// parameter does not exists
  template <class T>
  T getParameter(
      const std::string & name, const T & default_value,
      ParserParameterSearchCxt search_ctx = _ppsc_current_scope) const {
    try {
      T tmp = this->getParameter(name, search_ctx);
      return tmp;
    } catch (debug::Exception &) {
      return default_value;
    }
  }

  /* ------------------------------------------------------------------------ */
  /// Check if parameter exists within specified context
  bool hasParameter(
      const std::string & name,
      ParserParameterSearchCxt search_ctx = _ppsc_current_scope) const {

    Parameters::const_iterator it;
    if ((search_ctx & _ppsc_current_scope) != 0) {
      it = parameters.find(name);
    }

    if (it == parameters.end()) {
      if ((search_ctx & _ppsc_parent_scope) != 0 and
          parent_section != nullptr) {
        return parent_section->hasParameter(name, search_ctx);
      }
      return false;
    }
    return true;
  }

  /* --------------------------------------------------------------------------
   */
  /// Get value of given parameter in context
  template <class T>
  T getParameterValue(
      const std::string & name,
      ParserParameterSearchCxt search_ctx = _ppsc_current_scope) const {
    const ParserParameter & tmp_param = getParameter(name, search_ctx);
    T t = tmp_param;
    return t;
  }

  /* --------------------------------------------------------------------------
   */
  /// Get section name
  std::string getName() const { return name; }
  /// Get section type
  ParserType getType() const { return type; }
  /// Get section option
  std::string getOption(const std::string & def = "") const {
    return (not option.empty()) ? option : def;
  }

protected:
  void setParent(const ParserSection & sect) { parent_section = &sect; }

  /* ---------------------------------------------------------------------- */
  /* Members                                                                */
  /* ---------------------------------------------------------------------- */
private:
  /// Pointer to the parent section
  const ParserSection * parent_section{nullptr};
  /// Name of section
  std::string name;
  /// Type of section, see AKANTU_SECTION_TYPES
  ParserType type{ParserType::_not_defined};
  /// Section option
  std::string option;
  /// Map of parameters in section
  Parameters parameters;
  /// Multi-map of subsections
  SubSections sub_sections_by_type;
};

/* ------------------------------------------------------------------------ */
/* Parser Class                                                             */
/* ------------------------------------------------------------------------ */
/// Root of parsing tree, represents the global ParserSection
class Parser : public ParserSection {
public:
  Parser() : ParserSection("global", ParserType::_global) {}

  void parse(const std::string & filename);

  std::string getLastParsedFile() const;

  static bool isPermissive() { return permissive_parser; }

public:
  /// Parse real scalar
  static Real parseReal(const std::string & value,
                        const ParserSection & section);
  /// Parse real vector
  static Vector<Real> parseVector(const std::string & value,
                                  const ParserSection & section);
  /// Parse real matrix
  static Matrix<Real> parseMatrix(const std::string & value,
                                  const ParserSection & section);
  /// Parse real random parameter
  static RandomParameter<Real>
  parseRandomParameter(const std::string & value,
                       const ParserSection & section);

protected:
  /// General parse function
  template <class T, class Grammar>
  static T parseType(const std::string & value, Grammar & grammar);

protected:
  //  friend class Parsable;
  static bool permissive_parser;
  std::string last_parsed_file;
};

inline std::ostream & operator<<(std::ostream & stream,
                                 const ParserParameter & _this) {
  _this.printself(stream);
  return stream;
}

inline std::ostream & operator<<(std::ostream & stream,
                                 const ParserSection & section) {
  section.printself(stream);
  return stream;
}

} // namespace akantu

namespace std {
template <> struct iterator_traits<::akantu::Parser::const_section_iterator> {
  using iterator_category = input_iterator_tag;
  using value_type = ::akantu::ParserParameter;
  using difference_type = ptrdiff_t;
  using pointer = const ::akantu::ParserParameter *;
  using reference = const ::akantu::ParserParameter &;
};
} // namespace std

#include "parser_tmpl.hh"

#endif /* AKANTU_PARSER_HH_ */
