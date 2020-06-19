/**
 * @file   cppargparse_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Apr 03 2014
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  Implementation of the templated part of the commandline argument
 * parser
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
#include <sstream>
#include <stdexcept>

#ifndef __CPPARGPARSE_TMPL_HH__
#define __CPPARGPARSE_TMPL_HH__

namespace cppargparse {
/* -------------------------------------------------------------------------- */
/* Argument                                                                   */
/* -------------------------------------------------------------------------- */

/// internal description of arguments
struct ArgumentParser::_Argument : public Argument {
  _Argument() : Argument(), help(std::string()) {}
  ~_Argument() override = default;

  void setValues(std::vector<std::string> & values) {
    for (auto it = values.begin(); it != values.end(); ++it) {
      this->addValue(*it);
    }
  }

  virtual void addValue(std::string & value) = 0;
  virtual void setToDefault() = 0;
  virtual void setToConst() = 0;

  std::ostream & printDefault(std::ostream & stream) const {
    stream << std::boolalpha;
    if (has_default) {
      stream << " (default: ";
      this->_printDefault(stream);
      stream << ")";
    }
    if (has_const) {
      stream << " (const: ";
      this->_printConst(stream);
      stream << ")";
    }
    return stream;
  }

  virtual std::ostream & _printDefault(std::ostream & stream) const = 0;
  virtual std::ostream & _printConst(std::ostream & stream) const = 0;

  std::string help;
  int nargs{1};
  ArgumentType type{_string};
  bool required{false};
  bool parsed{false};
  bool has_default{false};
  bool has_const{false};
  std::vector<std::string> keys;
  bool is_positional{false};
};

/* -------------------------------------------------------------------------- */
/// typed storage of the arguments
template <class T>
class ArgumentParser::ArgumentStorage : public ArgumentParser::_Argument {
public:
  ArgumentStorage() : _default(T()), _const(T()), values(std::vector<T>()) {}

  void addValue(std::string & value) override {
    std::stringstream sstr(value);
    T t;
    sstr >> t;
    values.push_back(t);
  }

  void setToDefault() override {
    values.clear();
    values.push_back(_default);
  }

  void setToConst() override {
    values.clear();
    values.push_back(_const);
  }

  void printself(std::ostream & stream) const override {
    stream << this->name << " =";
    stream << std::boolalpha; // for boolean
    for (auto vit = this->values.begin(); vit != this->values.end(); ++vit) {
      stream << " " << *vit;
    }
  }

  std::ostream & _printDefault(std::ostream & stream) const override {
    stream << _default;
    return stream;
  }

  std::ostream & _printConst(std::ostream & stream) const override {
    stream << _const;
    return stream;
  }

  T _default;
  T _const;
  std::vector<T> values;
};

/* -------------------------------------------------------------------------- */
template <>
inline void
ArgumentParser::ArgumentStorage<std::string>::addValue(std::string & value) {
  values.push_back(value);
}

template <class T> struct is_vector {
  enum { value = false };
};

template <class T> struct is_vector<std::vector<T>> {
  enum { value = true };
};

/* -------------------------------------------------------------------------- */
template <class T, bool is_vector = cppargparse::is_vector<T>::value>
struct cast_helper {
  static T cast(const ArgumentParser::Argument & arg) {
    const auto & _arg =
        dynamic_cast<const ArgumentParser::ArgumentStorage<T> &>(arg);
    if (_arg.values.size() == 1) {
      return _arg.values[0];
    } else {
      throw std::length_error("Not enougth or too many argument where passed "
                              "for the command line argument: " +
                              arg.name);
    }
  }
};

template <class T> struct cast_helper<T, true> {
  static T cast(const ArgumentParser::Argument & arg) {
    const auto & _arg =
        dynamic_cast<const ArgumentParser::ArgumentStorage<T> &>(arg);
    return _arg.values;
  }
};

/* -------------------------------------------------------------------------- */
template <class T> ArgumentParser::Argument::operator T() const {
  return cast_helper<T>::cast(*this);
}

#if !defined(DOXYGEN)
template <> inline ArgumentParser::Argument::operator const std::string() const {
  return cast_helper<std::string>::cast(*this);
}

template <> inline ArgumentParser::Argument::operator unsigned int() const {
  return cast_helper<int>::cast(*this);
}
#endif

template <class T>
void ArgumentParser::addArgument(const std::string & name_or_flag,
                                 const std::string & help, int nargs,
                                 ArgumentType type, T def) {
  _Argument & arg = _addArgument(name_or_flag, help, nargs, type);
  dynamic_cast<ArgumentStorage<T> &>(arg)._default = def;
  arg.has_default = true;
}

template <class T>
void ArgumentParser::addArgument(const std::string & name_or_flag,
                                 const std::string & help, int nargs,
                                 ArgumentType type, T def, T cons) {
  _Argument & arg = _addArgument(name_or_flag, help, nargs, type);
  dynamic_cast<ArgumentStorage<T> &>(arg)._default = def;
  arg.has_default = true;
  dynamic_cast<ArgumentStorage<T> &>(arg)._const = cons;
  arg.has_const = true;
}

/* -------------------------------------------------------------------------- */
template <>
inline void
ArgumentParser::addArgument<const char *>(const std::string & name_or_flag,
                                          const std::string & help, int nargs,
                                          ArgumentType type, const char * def) {
  this->addArgument<std::string>(name_or_flag, help, nargs, type, def);
}

template <>
inline void
ArgumentParser::addArgument<unsigned int>(const std::string & name_or_flag,
                                          const std::string & help, int nargs,
                                          ArgumentType type, unsigned int def) {
  this->addArgument<int>(name_or_flag, help, nargs, type, def);
}

/* -------------------------------------------------------------------------- */
template <>
inline void ArgumentParser::addArgument<const char *>(
    const std::string & name_or_flag, const std::string & help, int nargs,
    ArgumentType type, const char * def, const char * cons) {
  this->addArgument<std::string>(name_or_flag, help, nargs, type, def, cons);
}

template <>
inline void ArgumentParser::addArgument<unsigned int>(
    const std::string & name_or_flag, const std::string & help, int nargs,
    ArgumentType type, unsigned int def, unsigned int cons) {
  this->addArgument<int>(name_or_flag, help, nargs, type, def, cons);
}
} // namespace cppargparse

#endif /* __AKANTU_CPPARGPARSE_TMPL_HH__ */
