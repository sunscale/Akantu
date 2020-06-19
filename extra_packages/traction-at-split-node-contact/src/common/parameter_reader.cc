/**
 * @file   parameter_reader.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Tue Dec 02 2014
 * @date last modification: Fri Feb 23 2018
 *
 * @brief  implementation of parameter reader
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
// std
#include <algorithm>
#include <fstream>
#include <iostream>
#include <utility>

// simtools
#include "parameter_reader.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
ParameterReader::ParameterReader()
    : data_types(), element_type_data(), string_data(), int_data(), uint_data(),
      real_data(), bool_data() {
  AKANTU_DEBUG_IN();

  // setup of types of element
  data_types.insert("elementtype");
  data_types.insert("string");
  data_types.insert("uint");
  data_types.insert("int");
  data_types.insert("real");
  data_types.insert("bool");
  //  data_types.insert("surface");

  // define conversion maps
  _input_to_akantu_element_types["_segment_2"] = akantu::_segment_2;
  _input_to_akantu_element_types["_segment_3"] = akantu::_segment_3;
  _input_to_akantu_element_types["_triangle_3"] = akantu::_triangle_3;
  _input_to_akantu_element_types["_triangle_6"] = akantu::_triangle_6;
  _input_to_akantu_element_types["_tetrahedron_4"] = akantu::_tetrahedron_4;
  _input_to_akantu_element_types["_tetrahedron_10"] = akantu::_tetrahedron_10;
  _input_to_akantu_element_types["_quadrangle_4"] = akantu::_quadrangle_4;
  _input_to_akantu_element_types["_quadrangle_8"] = akantu::_quadrangle_8;
  _input_to_akantu_element_types["_hexahedron_8"] = akantu::_hexahedron_8;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ParameterReader::readInputFile(std::string file_name) {
  AKANTU_DEBUG_IN();

  char comment_char = '#';
  char equal_char = '=';

  // open a file called file name
  std::ifstream infile;
  infile.open(file_name.c_str());

  if (!infile.good()) {
    std::cerr << "Cannot open file " << file_name << "!!!" << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string line;
  std::string clean_line;
  while (infile.good()) {
    getline(infile, line);
    clean_line = line;

    // take out comments
    size_t found_comment;
    found_comment = line.find_first_of(comment_char);
    if (found_comment != std::string::npos)
      clean_line = line.substr(0, found_comment);
    if (clean_line.empty())
      continue;

    std::stringstream sstr(clean_line);

    // check if data type exists
    std::string type;
    sstr >> type;
    std::transform(type.begin(), type.end(), type.begin(), ::tolower);
    if (this->data_types.find(type) == this->data_types.end()) {
      std::cerr << " *** WARNING *** Data type " << type << " does not exist"
                << " in this input data structure. Ignore line: ";
      std::cerr << clean_line << std::endl;
      continue;
    }

    std::string keyword;
    std::string equal;
    std::string value;

    // get keyword
    sstr >> keyword;
    size_t equal_p = keyword.find_first_of(equal_char);
    if (equal_p != std::string::npos) {
      equal = keyword.substr(equal_p, std::string::npos);
      keyword = keyword.substr(0, equal_p);
    }

    // get equal
    if (equal.empty())
      sstr >> equal;
    if (equal.length() != 1) {
      value = equal.substr(1, std::string::npos);
      equal = equal[0];
    }
    if (equal[0] != equal_char) {
      std::cerr << " *** WARNING *** Unrespected convention! Ignore line: ";
      std::cerr << clean_line << std::endl;
      continue;
    }

    // get value
    if (value.empty())
      sstr >> value;

    // no value
    if (value.empty()) {
      std::cerr << " *** WARNING *** No value given! Ignore line: ";
      std::cerr << clean_line << std::endl;
      continue;
    }

    // put value in map
    std::stringstream convert(value);
    if (type.compare("elementtype") == 0) {
      std::map<std::string, akantu::ElementType>::const_iterator it;
      it = this->_input_to_akantu_element_types.find(value);
      if (it != this->_input_to_akantu_element_types.end())
        this->element_type_data.insert(std::make_pair(keyword, it->second));
      else {
        std::cerr << " *** WARNING *** ElementType " << value
                  << " does not exist. Ignore line: ";
        std::cerr << clean_line << std::endl;
        continue;
      }
    } else if (type.compare("string") == 0) {
      this->string_data.insert(std::make_pair(keyword, value));
    }
    /*
    else if (type.compare("surface") == 0) {
      //Surface surf;
      UInt surf;
      convert >> surf;
      //this->surface_data.insert(std::make_pair(keyword,surf));
      this->uint_data.insert(std::make_pair(keyword,surf));
    }
    */
    else if (type.compare("int") == 0) {
      Int i;
      convert >> i;
      this->int_data.insert(std::make_pair(keyword, i));
    } else if (type.compare("uint") == 0) {
      UInt i;
      convert >> i;
      this->uint_data.insert(std::make_pair(keyword, i));
    } else if (type.compare("real") == 0) {
      Real r;
      convert >> r;
      this->real_data.insert(std::make_pair(keyword, r));
    } else if (type.compare("bool") == 0) {
      std::transform(value.begin(), value.end(), value.begin(), ::tolower);
      bool b;
      if (value.compare("true") == 0)
        b = true;
      else if (value.compare("false") == 0)
        b = false;
      else {
        std::cerr << " *** WARNING *** boolean cannot be " << value
                  << ". Ignore line: ";
        std::cerr << clean_line << std::endl;
        continue;
      }
      this->bool_data.insert(std::make_pair(keyword, b));
    } else {
      std::cerr << " *** ERROR *** Could not add data to InputData for line: ";
      std::cerr << clean_line << std::endl;
      continue;
      exit(EXIT_FAILURE);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ParameterReader::writeInputFile(std::string file_name) const {
  AKANTU_DEBUG_IN();

  // open file to write input information
  std::ofstream outfile;
  outfile.open(file_name.c_str());

  // element type
  for (std::map<std::string, akantu::ElementType>::const_iterator it =
           element_type_data.begin();
       it != element_type_data.end(); ++it) {
    for (std::map<std::string, ElementType>::const_iterator et =
             _input_to_akantu_element_types.begin();
         et != _input_to_akantu_element_types.end(); ++et) {
      if (it->second == et->second) {
        outfile << "ElementType " << it->first << " = " << et->first
                << std::endl;
        continue;
      }
    }
  }

  // string
  for (std::map<std::string, std::string>::const_iterator it =
           string_data.begin();
       it != string_data.end(); ++it)
    outfile << "string " << it->first << " = " << it->second << std::endl;

  // Surface
  /*
  for (std::map<std::string, akantu::Surface>::const_iterator it =
  surface_data.begin();
       it != surface_data.end(); ++it)
    outfile << "Surface " << it->first << " = " << it->second << std::endl;
  */
  // Int
  for (std::map<std::string, akantu::Int>::const_iterator it = int_data.begin();
       it != int_data.end(); ++it)
    outfile << "Int " << it->first << " = " << it->second << std::endl;

  // UInt
  for (std::map<std::string, akantu::UInt>::const_iterator it =
           uint_data.begin();
       it != uint_data.end(); ++it)
    outfile << "UInt " << it->first << " = " << it->second << std::endl;

  // Real
  for (std::map<std::string, akantu::Real>::const_iterator it =
           real_data.begin();
       it != real_data.end(); ++it)
    outfile << "Real " << it->first << " = " << it->second << std::endl;

  // Bool
  for (std::map<std::string, bool>::const_iterator it = bool_data.begin();
       it != bool_data.end(); ++it) {
    std::string b = "false";
    if (it->second)
      b = "true";
    outfile << "bool " << it->first << " = " << b << std::endl;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <>
akantu::UInt ParameterReader::get<akantu::UInt>(std::string key) const {
  std::map<std::string, akantu::UInt>::const_iterator it;
  it = this->uint_data.find(key);

  // if not in map
  if (it == this->uint_data.end()) {
    std::cerr << " *** ERROR *** This data was not in input file. "
              << "You need the following line in your input file: ";
    std::cerr << "UInt " << key << " = ???" << std::endl;
    exit(EXIT_FAILURE);
  }

  else
    return it->second;
}

/* -------------------------------------------------------------------------- */
template <>
akantu::ElementType
ParameterReader::get<akantu::ElementType>(std::string key) const {
  std::map<std::string, akantu::ElementType>::const_iterator it;
  it = this->element_type_data.find(key);

  // if not in map
  if (it == this->element_type_data.end()) {
    std::cerr << " *** ERROR *** This data was not in input file. "
              << "You need the following line in your input file: ";
    std::cerr << "ElementType " << key << " = ???" << std::endl;
    exit(EXIT_FAILURE);
  }

  else
    return it->second;
}

/* -------------------------------------------------------------------------- */
template <>
std::string ParameterReader::get<std::string>(std::string key) const {
  std::map<std::string, std::string>::const_iterator it;
  it = this->string_data.find(key);

  // if not in map
  if (it == this->string_data.end()) {
    std::cerr << " *** ERROR *** This data was not in input file. "
              << "You need the following line in your input file: ";
    std::cerr << "string " << key << " = ???" << std::endl;
    exit(EXIT_FAILURE);
  }

  else
    return it->second;
}

/* -------------------------------------------------------------------------- */
/*
template<>
akantu::Surface ParameterData::get<akantu::Surface>(std::string key) const {
  std::map<std::string,akantu::Surface>::const_iterator it;
  it = this->surface_data.find(key);

  // if not in map
  if (it == this->surface_data.end()) {
    std::cerr << " *** ERROR *** This data was not in input file. "
          << "You need the following line in your input file: ";
    std::cerr << "Surface " << key << " = ???" << std::endl;
    exit(EXIT_FAILURE);
  }

  else
    return it->second;
}
*/

/* -------------------------------------------------------------------------- */
template <>
akantu::Int ParameterReader::get<akantu::Int>(std::string key) const {
  std::map<std::string, akantu::Int>::const_iterator it;
  it = this->int_data.find(key);

  // if not in map
  if (it == this->int_data.end()) {
    std::cerr << " *** ERROR *** This data was not in input file. "
              << "You need the following line in your input file: ";
    std::cerr << "Int " << key << " = ???" << std::endl;
    exit(EXIT_FAILURE);
  }

  else
    return it->second;
}

/* -------------------------------------------------------------------------- */
template <>
akantu::Real ParameterReader::get<akantu::Real>(std::string key) const {
  std::map<std::string, akantu::Real>::const_iterator it;
  it = this->real_data.find(key);

  // if not in map
  if (it == this->real_data.end()) {
    std::cerr << " *** ERROR *** This data was not in input file. "
              << "You need the following line in your input file: ";
    std::cerr << "Real " << key << " = ???" << std::endl;
    exit(EXIT_FAILURE);
  }

  else
    return it->second;
}

/* -------------------------------------------------------------------------- */
template <> bool ParameterReader::get<bool>(std::string key) const {
  std::map<std::string, bool>::const_iterator it;
  it = this->bool_data.find(key);

  // if not in map
  if (it == this->bool_data.end()) {
    std::cerr << " *** ERROR *** This data was not in input file. "
              << "You need the following line in your input file: ";
    std::cerr << "bool " << key << " = ???" << std::endl;
    exit(EXIT_FAILURE);
  }

  else
    return it->second;
}

/* -------------------------------------------------------------------------- */
template <> bool ParameterReader::has<bool>(std::string key) const {
  std::map<std::string, bool>::const_iterator it;
  it = this->bool_data.find(key);
  return (it != this->bool_data.end());
}
template <> bool ParameterReader::has<std::string>(std::string key) const {
  std::map<std::string, std::string>::const_iterator it;
  it = this->string_data.find(key);
  return (it != this->string_data.end());
}
template <> bool ParameterReader::has<akantu::Int>(std::string key) const {
  std::map<std::string, akantu::Int>::const_iterator it;
  it = this->int_data.find(key);
  return (it != this->int_data.end());
}
template <> bool ParameterReader::has<akantu::UInt>(std::string key) const {
  std::map<std::string, akantu::UInt>::const_iterator it;
  it = this->uint_data.find(key);
  return (it != this->uint_data.end());
}
template <> bool ParameterReader::has<akantu::Real>(std::string key) const {
  std::map<std::string, akantu::Real>::const_iterator it;
  it = this->real_data.find(key);
  return (it != this->real_data.end());
}

/* -------------------------------------------------------------------------- */
void ParameterReader::printself(std::ostream & stream, int indent) const {
  AKANTU_DEBUG_IN();

  std::string space;
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
    ;

  stream << space << "ParameterReader [" << std::endl;
  /*
  stream << space << this->element_type_data << std::endl;
  stream << space << this->string_data << std::endl;
  stream << space << this->surface_data << std::endl;
  stream << space << this->int_data << std::endl;
  stream << space << this->uint_data << std::endl;
  stream << space << this->real_data << std::endl;
  stream << space << this->bool_data << std::endl;
  */
  stream << space << "]" << std::endl;

  AKANTU_DEBUG_OUT();
}

} // namespace akantu
