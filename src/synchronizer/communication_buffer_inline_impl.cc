/**
 * @file   communication_buffer_inline_impl.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Apr 14 2011
 * @date last modification: Thu Mar 27 2014
 *
 * @brief  CommunicationBuffer inline implementation
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
template<bool is_static>
inline void CommunicationBufferTemplated<is_static>::packResize(UInt size) {
  if(!is_static) {
    char * values = buffer.storage();
    buffer.resize(buffer.getSize() + size);
    ptr_pack = buffer.storage() + (ptr_pack - values);
    ptr_unpack = buffer.storage() + (ptr_unpack - values);
  }
}

/* -------------------------------------------------------------------------- */
template<bool is_static>
template<typename T>
inline CommunicationBufferTemplated<is_static> & CommunicationBufferTemplated<is_static>::operator<< (const T & to_pack) {
  packResize(sizeof(T));
  T * tmp = reinterpret_cast<T *>(ptr_pack);
  AKANTU_DEBUG_ASSERT(buffer.storage() + buffer.getSize() >= ptr_pack + sizeof(T),
		      "Packing too much data in the CommunicationBufferTemplated");
  *tmp = to_pack;
  ptr_pack += sizeof(T);
  return *this;
}

/* -------------------------------------------------------------------------- */
template<bool is_static>
template<typename T>
inline CommunicationBufferTemplated<is_static> & CommunicationBufferTemplated<is_static>::operator>> (T & to_unpack) {
  T * tmp = reinterpret_cast<T *>(ptr_unpack);
  to_unpack = *tmp;
  ptr_unpack += sizeof(T);
  return *this;
}


/* -------------------------------------------------------------------------- */
/* Specialization                                                             */
/* -------------------------------------------------------------------------- */

/**
 * Vector
 */

/* -------------------------------------------------------------------------- */
template<bool is_static>
template<typename T>
inline CommunicationBufferTemplated<is_static> & CommunicationBufferTemplated<is_static>::operator<< (const Vector<T> & to_pack) {
  UInt size = to_pack.size() * sizeof(T);
  packResize(size);
  AKANTU_DEBUG_ASSERT(buffer.storage() + buffer.getSize() >= ptr_pack + size,
		      "Packing too much data in the CommunicationBufferTemplated");
  memcpy(ptr_pack, to_pack.storage(), size);
  ptr_pack += size;
  return *this;
}

/* -------------------------------------------------------------------------- */
template<bool is_static>
template<typename T>
inline CommunicationBufferTemplated<is_static> & CommunicationBufferTemplated<is_static>::operator>> (Vector<T> & to_unpack) {
  UInt size = to_unpack.size() * sizeof(T);
  memcpy(to_unpack.storage(), ptr_unpack, size);
  ptr_unpack += size;
  return *this;
}

/**
 * Matrix
 */

/* -------------------------------------------------------------------------- */
template<bool is_static>
template<typename T>
inline CommunicationBufferTemplated<is_static> & CommunicationBufferTemplated<is_static>::operator<< (const Matrix<T> & to_pack) {
  UInt size = to_pack.size() * sizeof(Real);
  packResize(size);
  AKANTU_DEBUG_ASSERT(buffer.storage() + buffer.getSize() >= ptr_pack + size,
		      "Packing too much data in the CommunicationBufferTemplated");
  memcpy(ptr_pack, to_pack.storage(), size);
  ptr_pack += size;
  return *this;
}

/* -------------------------------------------------------------------------- */
template<bool is_static>
template<typename T>
inline CommunicationBufferTemplated<is_static> & CommunicationBufferTemplated<is_static>::operator>> (Matrix<T> & to_unpack) {
  UInt size = to_unpack.size() * sizeof(Real);
  memcpy(to_unpack.storage(), ptr_unpack, size);
  ptr_unpack += size;
  return *this;
}



/* -------------------------------------------------------------------------- */
template<bool is_static>
template<typename T>
inline void CommunicationBufferTemplated<is_static>::packIterable (T & to_pack) {
  UInt size = to_pack.size();
  operator<<(size);
  typename T::const_iterator it  = to_pack.begin();
  typename T::const_iterator end = to_pack.end();
  for(;it != end; ++it) operator<<(*it);
}

/* -------------------------------------------------------------------------- */
template<bool is_static>
template<typename T>
inline void CommunicationBufferTemplated<is_static>::unpackIterable (T & to_unpack) {
  UInt size;
  operator>>(size);
  to_unpack.resize(size);
  typename T::iterator it  = to_unpack.begin();
  typename T::iterator end = to_unpack.end();
  for(;it != end; ++it) operator>>(*it);
}

/**
 * std::vector<T>
 */
/* -------------------------------------------------------------------------- */

template<bool is_static>
template<typename T>
inline CommunicationBufferTemplated<is_static> & CommunicationBufferTemplated<is_static>::operator<< (const std::vector<T> & to_pack) {
  packIterable(to_pack);
  return *this;
}

/* -------------------------------------------------------------------------- */
template<bool is_static>
template<typename T>
inline CommunicationBufferTemplated<is_static> & CommunicationBufferTemplated<is_static>::operator>> (std::vector<T> & to_unpack) {
  unpackIterable(to_unpack);
  return *this;
}

/**
 * std::string
 */
/* -------------------------------------------------------------------------- */

template<bool is_static>
inline CommunicationBufferTemplated<is_static> & CommunicationBufferTemplated<is_static>::operator<< (const std::string & to_pack) {
  packIterable(to_pack);
  return *this;
}

/* -------------------------------------------------------------------------- */
template<bool is_static>
inline CommunicationBufferTemplated<is_static> & CommunicationBufferTemplated<is_static>::operator>> (std::string & to_unpack) {
  unpackIterable(to_unpack);
  return *this;
}

/* -------------------------------------------------------------------------- */
template<bool is_static>
template<typename T> inline std::string
CommunicationBufferTemplated<is_static>::extractStream(UInt block_size) {
  std::stringstream str;
  T * ptr = reinterpret_cast<T*>(buffer.storage());
  UInt sz = buffer.getSize()/sizeof(T);
  UInt sz_block = block_size/sizeof(T);

  UInt n_block = 0;
  for (UInt i = 0; i < sz; ++i) {
    if (i% sz_block == 0) {
      str << std::endl << n_block << " ";
      ++n_block;
    }
    str << *ptr << " ";
    ++ptr;
  }
  return str.str();
}

/* -------------------------------------------------------------------------- */
template<bool is_static>
inline void CommunicationBufferTemplated<is_static>::resize(UInt size) {
  if(!is_static) {
    buffer.resize(0);
  } else {
    buffer.resize(size);
  }
  reset();
#ifndef AKANTU_NDEBUG
  clear();
#endif
}

/* -------------------------------------------------------------------------- */
template<bool is_static>
inline void CommunicationBufferTemplated<is_static>::clear() {
  buffer.clear();
}

/* -------------------------------------------------------------------------- */
template<bool is_static>
inline void CommunicationBufferTemplated<is_static>::reset() {
  ptr_pack = buffer.storage();
  ptr_unpack = buffer.storage();
}

/* -------------------------------------------------------------------------- */
//template<bool is_static>
//inline CommunicationBufferTemplated<is_static> & CommunicationBufferTemplated<is_static>::packMeshData (const MeshData & to_pack, const ElementType & type) {

  //UInt size = to_pack.size();
  //operator<<(size);
  //typename std::vector<T>::iterator it  = to_pack.begin();
  //typename std::vector<T>::iterator end = to_pack.end();
  //for(;it != end; ++it) operator<<(*it);
  //return *this;

//}
