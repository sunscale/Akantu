/**
 * @file   communication_buffer.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  Buffer for packing and unpacking data
 *
 * @section LICENSE
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
#include "aka_array.hh"
#include "aka_common.hh"
#include "element.hh"

#ifndef __AKANTU_COMMUNICATION_BUFFER_HH__
#define __AKANTU_COMMUNICATION_BUFFER_HH__

namespace akantu {

template <bool is_static = true> class CommunicationBufferTemplated {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  explicit CommunicationBufferTemplated(UInt size) : buffer(size, 1, char()) {
    ptr_pack = buffer.storage();
    ptr_unpack = buffer.storage();
  };

  CommunicationBufferTemplated() : CommunicationBufferTemplated(0) {}

  CommunicationBufferTemplated(const CommunicationBufferTemplated & other) =
      delete;
  CommunicationBufferTemplated &
  operator=(const CommunicationBufferTemplated & other) = delete;

  CommunicationBufferTemplated(CommunicationBufferTemplated && other) = default;

  virtual ~CommunicationBufferTemplated() = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// reset to "empty"
  inline void reset();

  /// resize the internal buffer do not allocate on dynamic buffers
  inline void resize(UInt size);

  /// resize the internal buffer allocate always
  inline void reserve(UInt size);

  /// clear buffer context
  inline void clear();

private:
  inline void packResize(UInt size);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  inline char * storage() { return buffer.storage(); };
  inline const char * storage() const { return buffer.storage(); };
  /* ------------------------------------------------------------------------ */
  /* Operators                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// printing tool
  template <typename T> inline std::string extractStream(UInt packet_size);

  /// packing data
  template <typename T>
  inline CommunicationBufferTemplated & operator<<(const T & to_pack);

  template <typename T>
  inline CommunicationBufferTemplated & operator<<(const Vector<T> & to_pack);

  template <typename T>
  inline CommunicationBufferTemplated & operator<<(const Matrix<T> & to_pack);

  template <typename T>
  inline CommunicationBufferTemplated &
  operator<<(const std::vector<T> & to_pack);

  /// unpacking data
  template <typename T>
  inline CommunicationBufferTemplated & operator>>(T & to_unpack);

  template <typename T>
  inline CommunicationBufferTemplated & operator>>(Vector<T> & to_unpack);

  template <typename T>
  inline CommunicationBufferTemplated & operator>>(Matrix<T> & to_unpack);

  template <typename T>
  inline CommunicationBufferTemplated & operator>>(std::vector<T> & to_unpack);

  inline CommunicationBufferTemplated & operator<<(const std::string & to_pack);
  inline CommunicationBufferTemplated & operator>>(std::string & to_unpack);

private:
  template <typename T> inline void packIterable(T & to_pack);
  template <typename T> inline void unpackIterable(T & to_pack);

  /* ------------------------------------------------------------------------ */
  /* Accessor                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  template <typename T> static inline UInt sizeInBuffer(const T & data);
  template <typename T> static inline UInt sizeInBuffer(const Vector<T> & data);
  template <typename T> static inline UInt sizeInBuffer(const Matrix<T> & data);
  template <typename T>
  static inline UInt sizeInBuffer(const std::vector<T> & data);
  static inline UInt sizeInBuffer(const std::string & data);

  /// return the size in bytes of the stored values
  inline UInt getPackedSize() const { return ptr_pack - buffer.storage(); };
  /// return the size in bytes of data left to be unpacked
  inline UInt getLeftToUnpack() const {
    return buffer.size() - (ptr_unpack - buffer.storage());
  };
  /// return the global size allocated
  inline UInt size() const { return buffer.size(); };

  /// is the buffer empty
  inline bool empty() const {
    return (getPackedSize() == 0) and (getLeftToUnpack() == 0);
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// current position for packing
  char * ptr_pack;

  /// current position for unpacking
  char * ptr_unpack;

  /// storing buffer
  Array<char> buffer;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#if defined(AKANTU_INCLUDE_INLINE_IMPL)
#include "communication_buffer_inline_impl.cc"
#endif

using CommunicationBuffer = CommunicationBufferTemplated<true>;
using DynamicCommunicationBuffer = CommunicationBufferTemplated<false>;

} // namespace akantu

#endif /* __AKANTU_COMMUNICATION_BUFFER_HH__ */
