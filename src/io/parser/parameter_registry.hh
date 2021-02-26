/**
 * @file   parameter_registry.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Aug 09 2012
 * @date last modification: Tue Jan 30 2018
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
#include "parser.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_PARAMETER_REGISTRY_HH_
#define AKANTU_PARAMETER_REGISTRY_HH_

namespace akantu {
class ParserParameter;
}

namespace akantu {

/* -------------------------------------------------------------------------- */
/// Defines the access modes of parsable parameters
enum ParameterAccessType {
  _pat_internal = 0x0001,
  _pat_writable = 0x0010,
  _pat_readable = 0x0100,
  _pat_modifiable = 0x0110, //<_pat_readable | _pat_writable,
  _pat_parsable = 0x1000,
  _pat_parsmod = 0x1110 //< _pat_parsable | _pat_modifiable
};

/// Bit-wise operator between access modes
inline ParameterAccessType operator|(const ParameterAccessType & a,
                                     const ParameterAccessType & b) {
  auto tmp = ParameterAccessType(UInt(a) | UInt(b));
  return tmp;
}

/* -------------------------------------------------------------------------- */
template <typename T> class ParameterTyped;

/**
 * Interface for the Parameter
 */
class Parameter {
public:
  Parameter();
  Parameter(std::string name, std::string description,
            ParameterAccessType param_type);

  virtual ~Parameter() = default;
  /* ------------------------------------------------------------------------ */
  bool isInternal() const;
  bool isWritable() const;
  bool isReadable() const;
  bool isParsable() const;

  void setAccessType(ParameterAccessType ptype);

  /* ------------------------------------------------------------------------ */
  template <typename T, typename V> void set(const V & value);
  virtual void setAuto(const ParserParameter & param);
  template <typename T> T & get();
  template <typename T> const T & get() const;

  virtual inline operator Real() const { throw std::bad_cast(); };
  template <typename T> inline operator T() const;
  /* ------------------------------------------------------------------------ */
  virtual void printself(std::ostream & stream) const;

  virtual const std::type_info & type() const = 0;

protected:
  /// Returns const instance of templated sub-class ParameterTyped
  template <typename T> const ParameterTyped<T> & getParameterTyped() const;

  /// Returns instance of templated sub-class ParameterTyped
  template <typename T> ParameterTyped<T> & getParameterTyped();

protected:
  /// Name of parameter
  std::string name;

private:
  /// Description of parameter
  std::string description;
  /// Type of access
  ParameterAccessType param_type{_pat_internal};
};

/* -------------------------------------------------------------------------- */
/* Typed Parameter                                                            */
/* -------------------------------------------------------------------------- */
/**
 * Type parameter transfering a ParserParameter (string: string) to a typed
 * parameter in the memory of the p
 */
template <typename T> class ParameterTyped : public Parameter {
public:
  ParameterTyped(const std::string & name, const std::string & description,
                 ParameterAccessType param_type, T & param);

  /* ------------------------------------------------------------------------ */
  template <typename V> void setTyped(const V & value);
  void setAuto(const ParserParameter & value) override;
  T & getTyped();
  const T & getTyped() const;

  void printself(std::ostream & stream) const override;

  inline operator Real() const override;

  inline const std::type_info & type() const override { return typeid(T); }

private:
  /// Value of parameter
  T & param;
};

/* -------------------------------------------------------------------------- */
/* Parsable Interface                                                         */
/* -------------------------------------------------------------------------- */
/// Defines interface for classes to manipulate parsable parameters
class ParameterRegistry {
public:
  ParameterRegistry();
  virtual ~ParameterRegistry();

  /* ------------------------------------------------------------------------ */
  /// Add parameter to the params map
  template <typename T>
  void registerParam(const std::string & name, T & variable, ParameterAccessType type,
                     const std::string & description = "");

  /// Add parameter to the params map (with default value)
  template <typename T>
  void registerParam(const std::string &name, T & variable, const T & default_value,
                     ParameterAccessType type,
                     const std::string & description = "");

  /*------------------------------------------------------------------------- */
protected:
  void registerSubRegistry(const ID & id, ParameterRegistry & registry);

  /* ------------------------------------------------------------------------ */
public:
  /// Set value to a parameter (with possible different type)
  template <typename T, typename V>
  void setMixed(const std::string & name, const V & value);

  /// Set value to a parameter
  template <typename T> void set(const std::string & name, const T & value);

  /// Get value of a parameter
  inline const Parameter & get(const std::string & name) const;

  /// Get value of a parameter
  inline Parameter & get(const std::string & name);

  std::vector<ID> listParameters() const {
    std::vector<ID> params;
    for (const auto & pair : this->params) {
      params.push_back(pair.first);
    }
    return params;
  }

  std::vector<ID> listSubRegisteries() const {
    std::vector<ID> subs;
    for (const auto & pair : this->sub_registries) {
      subs.push_back(pair.first);
    }
    return subs;
  }

protected:
  template <typename T> T & get_(const std::string & name);

protected:
  void setParameterAccessType(const std::string & name,
                              ParameterAccessType ptype);
  /* ------------------------------------------------------------------------ */
  virtual void printself(std::ostream & stream, int indent) const;

protected:
  /// Parameters map
  using Parameters = std::map<std::string, Parameter *>;
  Parameters params;

  /// list of sub-registries
  using SubRegisteries = std::map<std::string, ParameterRegistry *>;
  SubRegisteries sub_registries;

  /// should accessor check in sub registries
  bool consisder_sub{true};
};

} // namespace akantu

#include "parameter_registry_tmpl.hh"

#endif /* AKANTU_PARAMETER_REGISTRY_HH_ */
