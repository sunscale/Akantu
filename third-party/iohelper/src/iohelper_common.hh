/**
 * @file   iohelper_common.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Mar 11 2010
 * @date last modification: Thu Oct 10 2013
 *
 * @brief  header for common types
 *
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * IOHelper is free  software: you can redistribute it and/or  modify it under
 * the terms  of the  GNU Lesser  General Public  License as  published by  the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * IOHelper is  distributed in the  hope that it  will be useful, but  WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for
 * more details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with IOHelper. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#ifndef IOHELPER_COMMON_H_
#define IOHELPER_COMMON_H_
/* -------------------------------------------------------------------------- */

#define USING_ZLIB
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <string>

/* -------------------------------------------------------------------------- */

namespace iohelper {

using UInt = unsigned int;
using Int = int;
using Real = double;

/* -------------------------------------------------------------------------- */
enum DataType { _bool, _uint, _int, _float, _double, _int64, _uint64, _uint8 };

enum IndexingMode { C_MODE = 0, FORTRAN_MODE = 1 };

/* -------------------------------------------------------------------------- */

#if __cplusplus <= 199711L
enum ElemType {
#else
enum ElemType : unsigned int {
#endif
  TRIANGLE1,
  TRIANGLE2,
  TETRA1,
  TETRA2,
  POINT_SET,
  LINE1,
  LINE2,
  QUAD1,
  QUAD2,
  HEX1,
  HEX2,
  BEAM2,
  BEAM3,
  PRISM1,
  PRISM2,
  COH1D2,
  COH2D4,
  COH2D6,
  COH3D6,
  COH3D12,
  COH3D8,
  MAX_ELEM_TYPE
};

/* -------------------------------------------------------------------------- */

enum FileStorageMode { TEXT = 0, BASE64 = 1, COMPRESSED = 2 };

enum TextDumpMode { _tdm_space, _tdm_csv };

/* -------------------------------------------------------------------------- */
static UInt nb_node_per_elem[MAX_ELEM_TYPE] __attribute__((unused)) = {
    3,  // TRIANGLE1
    6,  // TRIANGLE2
    4,  // TETRA1
    10, // TETRA2
    1,  // POINT_SET
    2,  // LINE1
    3,  // LINE2
    4,  // QUAD1
    8,  // QUAD2
    8,  // HEX1
    20, // HEX2
    2,  // BEAM2
    2,  // BEAM3
    6,  // PRISM1
    15, // PRISM2
    2,  // COH1D2
    4,  // COH2D4
    6,  // COH2D6
    6,  // COH3D6
    12, // COH3D12
    8,  // COH3D8
};

/* -------------------------------------------------------------------------- */
static UInt nb_quad_points[MAX_ELEM_TYPE] __attribute__((unused)) = {
    1,  // TRIANGLE1
    3,  // TRIANGLE2
    1,  // TETRA1
    4,  // TETRA2
    0,  // POINT_SET
    1,  // LINE1
    2,  // LINE2
    4,  // QUAD1
    9,  // QUAD2
    8,  // HEX1
    27, // HEX2
    2,  // BEAM2
    3,  // BEAM3
    6,  // PRISM1
    8,  // PRISM2
    1,  // COH1D2
    1,  // COH2D4
    2,  // COH2D6
    1,  // COH3D6
    3,  // COH3D12
    1,  // COH3D8
};

/* -------------------------------------------------------------------------- */
template <typename T> class IOHelperVector {
public:
  virtual ~IOHelperVector() = default;

  inline IOHelperVector(T * ptr, UInt size) {
    this->ptr = ptr;
    this->_size = size;
  };

  inline UInt size() const { return _size; };

  inline const T & operator[](UInt i) const { return ptr[i]; };

  inline const T * getPtr() const { return ptr; };

private:
  T * ptr;
  UInt _size;
};

/* -------------------------------------------------------------------------- */
/* Iterator interface                                                         */
/* -------------------------------------------------------------------------- */
template <typename T, class daughter, class ret_cont = IOHelperVector<T>>
class iterator {
public:
  using type = ret_cont;

  virtual ~iterator() = default;
  virtual bool operator!=(const daughter & it) const = 0;
  virtual daughter & operator++() = 0;
  virtual ret_cont operator*() = 0;
  //! This function is only for the element iterators
  virtual ElemType element_type() { return MAX_ELEM_TYPE; }
};

/* -------------------------------------------------------------------------- */

class IOHelperException : public std::exception {

public:
  enum ErrorType {
    _et_non_homogeneous_data,
    _et_unknown_visitor_stage,
    _et_file_error,
    _et_missing_field,
    _et_data_type,
    _et_options_error
  };

public:
  IOHelperException(const std::string & message,
                    const ErrorType type) noexcept {
    this->message = message;
    this->type = type;
  };

  ~IOHelperException() noexcept override = default;

  const char * what() const noexcept override { return message.c_str(); };

private:
  std::string message;
  ErrorType type;
};

/* -------------------------------------------------------------------------- */
#define IOHELPER_THROW(x, type)                                                \
  {                                                                            \
    std::stringstream ioh_throw_sstr;                                          \
    ioh_throw_sstr << __FILE__ << ":" << __LINE__ << ":"                       \
                   << __PRETTY_FUNCTION__ << ": " << x; /* NOLINT */           \
    std::string ioh_message(ioh_throw_sstr.str());                             \
    throw ::iohelper::IOHelperException(ioh_message,                           \
                                        ::iohelper::IOHelperException::type);  \
  }

/* -------------------------------------------------------------------------- */
template <typename T> DataType getDataType();

#define DEFINE_GET_DATA_TYPE(type, data_type)                                  \
  template <> inline DataType getDataType<type>() { return data_type; }

DEFINE_GET_DATA_TYPE(bool, _bool)
DEFINE_GET_DATA_TYPE(ElemType, _int)
DEFINE_GET_DATA_TYPE(int, _int)
DEFINE_GET_DATA_TYPE(unsigned int, _uint)
DEFINE_GET_DATA_TYPE(float, _float)
DEFINE_GET_DATA_TYPE(double, _double)
DEFINE_GET_DATA_TYPE(long int, _int64)
DEFINE_GET_DATA_TYPE(unsigned long int, _uint64)
DEFINE_GET_DATA_TYPE(std::uint8_t, _uint8)

#undef DEFINE_GET_DATA_TYPE

inline std::ostream & operator<<(std::ostream & stream, DataType type) {
  switch (type) {
  case _bool:
    stream << "bool";
    break;
  case _uint:
    stream << "uint32";
    break;
  case _int:
    stream << "int32";
    break;
  case _float:
    stream << "float32";
    break;
  case _double:
    stream << "float64";
    break;
  case _uint64:
    stream << "uint64";
    break;
  case _int64:
    stream << "int64";
    break;
  case _uint8:
    stream << "uint8";
    break;
  }
  return stream;
}

inline std::ostream & operator<<(std::ostream & stream, TextDumpMode mode) {
  switch (mode) {
  case _tdm_space:
    stream << "space";
    break;
  case _tdm_csv:
    stream << "csv";
    break;
  }
  return stream;
}

} // namespace iohelper

#endif /* IOHELPER_COMMON_H_ */
