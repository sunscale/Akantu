/**
 * @file   parameter_reader.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 *
 * @brief  for simulations to read parameters from an input file
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#ifndef __AST_PARAMETER_READER_HH__
#define __AST_PARAMETER_READER_HH__

/* -------------------------------------------------------------------------- */
// std
#include <map>
#include <set>

// akantu
#include "aka_common.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
class ParameterReader {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  ParameterReader();
  virtual ~ParameterReader(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// read input file
  void readInputFile(std::string file_name);

  /// write input file
  void writeInputFile(std::string file_name) const;

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  ///
  template <typename T> T get(std::string key) const;

  template <typename T> bool has(std::string key) const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// type of data available
  std::set<std::string> data_types;

  /// data
  std::map<std::string, akantu::ElementType> element_type_data;
  std::map<std::string, std::string> string_data;
  std::map<std::string, akantu::Int> int_data;
  std::map<std::string, akantu::UInt> uint_data;
  std::map<std::string, akantu::Real> real_data;
  std::map<std::string, bool> bool_data;

  /// convert string to element type
  std::map<std::string, ElementType> _input_to_akantu_element_types;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "parameter_reader_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream,
                                 const ParameterReader & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#endif /* __AST_PARAMETER_READER_HH__ */
