/**
 * @file   communication_buffer_inline_impl.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Apr 14 2011
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  CommunicationBuffer inline implementation
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
template <bool is_static>
template <typename T>
inline UInt CommunicationBufferTemplated<is_static>::sizeInBuffer(const T &) {
  return sizeof(T);
}

template <bool is_static>
template <typename T>
inline UInt
CommunicationBufferTemplated<is_static>::sizeInBuffer(const Vector<T> & data) {
  UInt size = data.size() * sizeof(T);
  return size;
}

template <bool is_static>
template <typename T>
inline UInt
CommunicationBufferTemplated<is_static>::sizeInBuffer(const Matrix<T> & data) {
  UInt size = data.size() * sizeof(T);
  return size;
}

template <bool is_static>
template <typename T>
inline UInt CommunicationBufferTemplated<is_static>::sizeInBuffer(
    const std::vector<T> & data) {
  UInt size = data.size() * sizeof(T) + sizeof(size_t);
  return size;
}

template <bool is_static>
inline UInt CommunicationBufferTemplated<is_static>::sizeInBuffer(
    const std::string & data) {
  UInt size = data.size() * sizeof(std::string::value_type) + sizeof(size_t);
  return size;
}

/* -------------------------------------------------------------------------- */
template <bool is_static>
inline void CommunicationBufferTemplated<is_static>::packResize(UInt size) {
  if (not is_static) {
    char * values = buffer.storage();
    auto nb_packed = ptr_pack - values;

    if (buffer.size() > nb_packed + size)
      return;

    buffer.resize(nb_packed + size);
    ptr_pack = buffer.storage() + nb_packed;
    ptr_unpack = buffer.storage() + (ptr_unpack - values);
  }
}

/* -------------------------------------------------------------------------- */
template <bool is_static>
template <typename T>
inline CommunicationBufferTemplated<is_static> &
CommunicationBufferTemplated<is_static>::operator<<(const T & to_pack) {
  UInt size = sizeInBuffer(to_pack);
  packResize(size);
  AKANTU_DEBUG_ASSERT(
      (buffer.storage() + buffer.size()) >= (ptr_pack + size),
      "Packing too much data in the CommunicationBufferTemplated");
  memcpy(ptr_pack, reinterpret_cast<const char *>(&to_pack), size);
  ptr_pack += size;
  return *this;
}

/* -------------------------------------------------------------------------- */
template <bool is_static>
template <typename T>
inline CommunicationBufferTemplated<is_static> &
CommunicationBufferTemplated<is_static>::operator>>(T & to_unpack) {
  UInt size = sizeInBuffer(to_unpack);

  alignas(alignof(T)) std::array<char, sizeof(T)> aligned_ptr;
  memcpy(aligned_ptr.data(), ptr_unpack, size);

  auto * tmp = reinterpret_cast<T *>(aligned_ptr.data());
  AKANTU_DEBUG_ASSERT(
      (buffer.storage() + buffer.size()) >= (ptr_unpack + size),
      "Unpacking too much data in the CommunicationBufferTemplated");
  to_unpack = *tmp;
  // memcpy(reinterpret_cast<char *>(&to_unpack), ptr_unpack, size);
  ptr_unpack += size;
  return *this;
}

/* -------------------------------------------------------------------------- */
/* Specialization                                                             */
/* -------------------------------------------------------------------------- */

/**
 * Vector
 */

/* -------------------------------------------------------------------------- */
template <bool is_static>
template <typename T>
inline CommunicationBufferTemplated<is_static> &
CommunicationBufferTemplated<is_static>::operator<<(const Vector<T> & to_pack) {
  UInt size = sizeInBuffer(to_pack);
  packResize(size);
  AKANTU_DEBUG_ASSERT(
      (buffer.storage() + buffer.size()) >= (ptr_pack + size),
      "Packing too much data in the CommunicationBufferTemplated");
  memcpy(ptr_pack, to_pack.storage(), size);
  ptr_pack += size;
  return *this;
}

/* -------------------------------------------------------------------------- */
template <bool is_static>
template <typename T>
inline CommunicationBufferTemplated<is_static> &
CommunicationBufferTemplated<is_static>::operator>>(Vector<T> & to_unpack) {
  UInt size = sizeInBuffer(to_unpack);
  AKANTU_DEBUG_ASSERT(
      (buffer.storage() + buffer.size()) >= (ptr_unpack + size),
      "Unpacking too much data in the CommunicationBufferTemplated");
  memcpy(to_unpack.storage(), ptr_unpack, size);
  ptr_unpack += size;
  return *this;
}

/**
 * Matrix
 */

/* -------------------------------------------------------------------------- */
template <bool is_static>
template <typename T>
inline CommunicationBufferTemplated<is_static> &
CommunicationBufferTemplated<is_static>::operator<<(const Matrix<T> & to_pack) {
  UInt size = sizeInBuffer(to_pack);
  packResize(size);
  AKANTU_DEBUG_ASSERT(
      (buffer.storage() + buffer.size()) >= (ptr_pack + size),
      "Packing too much data in the CommunicationBufferTemplated");
  memcpy(ptr_pack, to_pack.storage(), size);
  ptr_pack += size;
  return *this;
}

/* -------------------------------------------------------------------------- */
template <bool is_static>
template <typename T>
inline CommunicationBufferTemplated<is_static> &
CommunicationBufferTemplated<is_static>::operator>>(Matrix<T> & to_unpack) {
  UInt size = sizeInBuffer(to_unpack);
  AKANTU_DEBUG_ASSERT(
      (buffer.storage() + buffer.size()) >= (ptr_unpack + size),
      "Unpacking too much data in the CommunicationBufferTemplated");
  memcpy(to_unpack.storage(), ptr_unpack, size);
  ptr_unpack += size;
  return *this;
}

/* -------------------------------------------------------------------------- */
template <bool is_static>
template <typename T>
inline void CommunicationBufferTemplated<is_static>::packIterable(T & to_pack) {
  operator<<(size_t(to_pack.size()));
  auto it = to_pack.begin();
  auto end = to_pack.end();
  for (; it != end; ++it)
    operator<<(*it);
}

/* -------------------------------------------------------------------------- */
template <bool is_static>
template <typename T>
inline void
CommunicationBufferTemplated<is_static>::unpackIterable(T & to_unpack) {
  size_t size;
  operator>>(size);
  to_unpack.resize(size);
  auto it = to_unpack.begin();
  auto end = to_unpack.end();
  for (; it != end; ++it)
    operator>>(*it);
}

/**
 * std::vector<T>
 */
/* -------------------------------------------------------------------------- */

template <bool is_static>
template <typename T>
inline CommunicationBufferTemplated<is_static> &
CommunicationBufferTemplated<is_static>::
operator<<(const std::vector<T> & to_pack) {
  packIterable(to_pack);
  return *this;
}

/* -------------------------------------------------------------------------- */
template <bool is_static>
template <typename T>
inline CommunicationBufferTemplated<is_static> &
CommunicationBufferTemplated<is_static>::
operator>>(std::vector<T> & to_unpack) {
  unpackIterable(to_unpack);
  return *this;
}

/**
 * std::string
 */
/* -------------------------------------------------------------------------- */

template <bool is_static>
inline CommunicationBufferTemplated<is_static> &
CommunicationBufferTemplated<is_static>::
operator<<(const std::string & to_pack) {
  packIterable(to_pack);
  return *this;
}

/* -------------------------------------------------------------------------- */
template <bool is_static>
inline CommunicationBufferTemplated<is_static> &
CommunicationBufferTemplated<is_static>::operator>>(std::string & to_unpack) {
  unpackIterable(to_unpack);
  return *this;
}

/* -------------------------------------------------------------------------- */
template <bool is_static>
template <typename T>
inline std::string
CommunicationBufferTemplated<is_static>::extractStream(UInt block_size) {
  std::stringstream str;
  auto * ptr = reinterpret_cast<T *>(buffer.storage());
  UInt sz = buffer.size() / sizeof(T);
  UInt sz_block = block_size / sizeof(T);

  UInt n_block = 0;
  for (UInt i = 0; i < sz; ++i) {
    if (i % sz_block == 0) {
      str << std::endl << n_block << " ";
      ++n_block;
    }
    str << *ptr << " ";
    ++ptr;
  }
  return str.str();
}

/* -------------------------------------------------------------------------- */
template <bool is_static>
inline void CommunicationBufferTemplated<is_static>::resize(UInt size) {
  if (!is_static) {
    buffer.resize(0, 0);
  } else {
    buffer.resize(size, 0);
  }
  reset();
#ifndef AKANTU_NDEBUG
  clear();
#endif
}

/* -------------------------------------------------------------------------- */
template <bool is_static>
inline void CommunicationBufferTemplated<is_static>::reserve(UInt size) {
  char * values = buffer.storage();
  auto nb_packed = ptr_pack - values;

  buffer.resize(size);
  ptr_pack = buffer.storage() + nb_packed;
  ptr_unpack = buffer.storage() + (ptr_unpack - values);
}

/* -------------------------------------------------------------------------- */
template <bool is_static>
inline void CommunicationBufferTemplated<is_static>::clear() {
  buffer.clear();
}

/* -------------------------------------------------------------------------- */
template <bool is_static>
inline void CommunicationBufferTemplated<is_static>::reset() {
  ptr_pack = buffer.storage();
  ptr_unpack = buffer.storage();
}

/* -------------------------------------------------------------------------- */
// template<bool is_static>
// inline CommunicationBufferTemplated<is_static> &
// CommunicationBufferTemplated<is_static>::packMeshData (const MeshData &
// to_pack, const ElementType & type) {

// UInt size = to_pack.size();
// operator<<(size);
// typename std::vector<T>::iterator it  = to_pack.begin();
// typename std::vector<T>::iterator end = to_pack.end();
// for(;it != end; ++it) operator<<(*it);
// return *this;

//}
