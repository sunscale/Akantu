/**
 * @file   parameter_reader.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Fri Feb 23 2018
 *
 * @brief  for simulations to read parameters from an input file
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#ifndef AST_PARAMETER_READER_HH_
#define AST_PARAMETER_READER_HH_

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

//#include "parameter_reader_inline_impl.hh"

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream,
                                 const ParameterReader & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#endif /* AST_PARAMETER_READER_HH_ */
