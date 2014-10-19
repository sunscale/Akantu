/**
 * @file   element_type_map_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Aug 31 2011
 * @date last modification: Thu Jun 05 2014
 *
 * @brief  implementation of template functions of the ElementTypeMap and
 * ElementTypeMapArray classes
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef __AKANTU_ELEMENT_TYPE_MAP_TMPL_HH__
#define __AKANTU_ELEMENT_TYPE_MAP_TMPL_HH__

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/* ElementTypeMap                                                              */
/* -------------------------------------------------------------------------- */
template<class Stored, typename SupportType>
inline std::string ElementTypeMap<Stored, SupportType>::printType(const SupportType & type,
                                                                 const GhostType & ghost_type) {
  std::stringstream sstr; sstr << "(" << ghost_type << ":" << type << ")";
  return sstr.str();
}

/* -------------------------------------------------------------------------- */
template<class Stored, typename SupportType>
inline bool ElementTypeMap<Stored, SupportType>::exists(const SupportType & type, const GhostType & ghost_type) const {
  return this->getData(ghost_type).find(type) != this->getData(ghost_type).end();
}

/* -------------------------------------------------------------------------- */
template<class Stored, typename SupportType>
inline const Stored & ElementTypeMap<Stored, SupportType>::operator()(const SupportType & type,
                                                                     const GhostType & ghost_type) const {
  typename DataMap::const_iterator it =
    this->getData(ghost_type).find(type);

  if(it == this->getData(ghost_type).end())
    AKANTU_EXCEPTION("No element of type "
                     << ElementTypeMap::printType(type, ghost_type)
                     << " in this ElementTypeMap<"
                     << debug::demangle(typeid(Stored).name()) << "> class");
  return it->second;
}

/* -------------------------------------------------------------------------- */
template<class Stored, typename SupportType>
inline Stored & ElementTypeMap<Stored, SupportType>::operator()(const SupportType & type,
							       const GhostType & ghost_type) {
  return this->getData(ghost_type)[type];
}

/* -------------------------------------------------------------------------- */
template<class Stored, typename SupportType>
inline Stored & ElementTypeMap<Stored, SupportType>::operator()(const Stored & insert,
                                                               const SupportType & type,
                                                               const GhostType & ghost_type) {
  typename DataMap::iterator it =
    this->getData(ghost_type).find(type);

  if(it != this->getData(ghost_type).end()) {
    AKANTU_EXCEPTION("Element of type "
                     << ElementTypeMap::printType(type, ghost_type)
                     << " already in this ElementTypeMap<"
                     << debug::demangle(typeid(Stored).name()) << "> class");
  } else {
    DataMap & data = this->getData(ghost_type);
    const std::pair<typename DataMap::iterator, bool> & res =
      data.insert(std::pair<ElementType, Stored>(type, insert));
    it = res.first;
  }

  return it->second;
}


/* -------------------------------------------------------------------------- */
template<class Stored, typename SupportType>
inline typename ElementTypeMap<Stored, SupportType>::DataMap &
ElementTypeMap<Stored, SupportType>::getData(GhostType ghost_type) {
  if(ghost_type == _not_ghost) return data;
  else return ghost_data;
}

/* -------------------------------------------------------------------------- */
template<class Stored, typename SupportType>
inline const typename ElementTypeMap<Stored, SupportType>::DataMap &
ElementTypeMap<Stored, SupportType>::getData(GhostType ghost_type) const {
  if(ghost_type == _not_ghost) return data;
  else return ghost_data;
}

/* -------------------------------------------------------------------------- */
/// Works only if stored is a pointer to a class with a printself method
template<class Stored, typename SupportType>
void ElementTypeMap<Stored, SupportType>::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "ElementTypeMap<" << debug::demangle(typeid(Stored).name()) << "> [" << std::endl;
  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;

    const DataMap & data = getData(gt);
    typename DataMap::const_iterator it;
    for(it = data.begin(); it != data.end(); ++it) {
      stream << space << space << ElementTypeMap::printType(it->first, gt) << std::endl;
    }
  }
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
template<class Stored, typename SupportType>
ElementTypeMap<Stored, SupportType>::ElementTypeMap() {
  AKANTU_DEBUG_IN();

  // std::stringstream sstr;
  // if(parent_id != "") sstr << parent_id << ":";
  // sstr << id;

  // this->id = sstr.str();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<class Stored, typename SupportType>
ElementTypeMap<Stored, SupportType>::~ElementTypeMap() {

}

/* -------------------------------------------------------------------------- */
/* ElementTypeMapArray                                                        */
/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType>
inline Array<T> & ElementTypeMapArray<T, SupportType>::alloc(UInt size,
                                                            UInt nb_component,
                                                            const SupportType & type,
                                                            const GhostType & ghost_type) {
  std::string ghost_id = "";
  if (ghost_type == _ghost) ghost_id = ":ghost";

  Array<T> * tmp;

  typename ElementTypeMapArray<T, SupportType>::DataMap::iterator it =
    this->getData(ghost_type).find(type);

  if(it == this->getData(ghost_type).end()) {
    std::stringstream sstr; sstr << this->id << ":" << type << ghost_id;
    tmp = &(Memory::alloc<T>(sstr.str(), size,
			     nb_component, T()));
    std::stringstream sstrg; sstrg << ghost_type;
    //tmp->setTag(sstrg.str());
    this->getData(ghost_type)[type] = tmp;
  } else {
    AKANTU_DEBUG_INFO("The vector " << this->id << this->printType(type, ghost_type)
		      << " already exists, it is resized instead of allocated.");
    tmp = it->second;
    it->second->resize(size);
  }

  return *tmp;
}

/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType>
inline void ElementTypeMapArray<T, SupportType>::alloc(UInt size,
                                                      UInt nb_component,
                                                      const SupportType & type) {
  this->alloc(size, nb_component, type, _not_ghost);
  this->alloc(size, nb_component, type, _ghost);
}

/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType>
inline void ElementTypeMapArray<T, SupportType>::free() {
  AKANTU_DEBUG_IN();

  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;

    DataMap & data = this->getData(gt);
    typename DataMap::const_iterator it;
    for(it = data.begin(); it != data.end(); ++it) {
      dealloc(it->second->getID());
    }
    data.clear();
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType>
inline const Array<T> & ElementTypeMapArray<T, SupportType>::operator()(const SupportType & type,
                                                                       const GhostType & ghost_type) const {
  typename ElementTypeMapArray<T, SupportType>::DataMap::const_iterator it =
    this->getData(ghost_type).find(type);

  if(it == this->getData(ghost_type).end())
    AKANTU_EXCEPTION("No element of type "
                     << ElementTypeMapArray::printType(type, ghost_type)
                     << " in this const ElementTypeMapArray<"
                     << debug::demangle(typeid(T).name()) << "> class(\""
                     << this->id << "\")");
  return *(it->second);
}

/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType>
inline Array<T> & ElementTypeMapArray<T, SupportType>::operator()(const SupportType & type,
                                                                 const GhostType & ghost_type) {
  typename ElementTypeMapArray<T, SupportType>::DataMap::iterator it =
    this->getData(ghost_type).find(type);

  if(it == this->getData(ghost_type).end())
    AKANTU_EXCEPTION("No element of type "
                     << ElementTypeMapArray::printType(type, ghost_type)
                     << " in this ElementTypeMapArray<"
                     << debug::demangle(typeid(T).name()) << "> class (\""
                     << this->id << "\")");

  return *(it->second);
}

/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType>
inline void ElementTypeMapArray<T, SupportType>::setArray(const SupportType & type,
                                                         const GhostType & ghost_type,
                                                         const Array<T> & vect) {
  typename ElementTypeMapArray<T, SupportType>::DataMap::iterator it =
    this->getData(ghost_type).find(type);

  if(AKANTU_DEBUG_TEST(dblWarning) && it != this->getData(ghost_type).end() && it->second != &vect) {
    AKANTU_DEBUG_WARNING("The Array " << this->printType(type, ghost_type)
                         << " is already registred, this call can lead to a memory leak.");
  }

  this->getData(ghost_type)[type] = &(const_cast<Array<T> &>(vect));
}

/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType>
inline void ElementTypeMapArray<T, SupportType>::onElementsRemoved(const ElementTypeMapArray<UInt> & new_numbering) {
  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;
    ElementTypeMapArray<UInt>::type_iterator it  = new_numbering.firstType(_all_dimensions, gt, _ek_not_defined);
    ElementTypeMapArray<UInt>::type_iterator end = new_numbering.lastType(_all_dimensions, gt, _ek_not_defined);
    for (; it != end; ++it) {
      SupportType type = *it;
      if(this->exists(type, gt)){
        const Array<UInt> & renumbering = new_numbering(type, gt);
        Array<T> & vect = this->operator()(type, gt);
        UInt nb_component = vect.getNbComponent();
        Array<T> tmp(renumbering.getSize(), nb_component);
        UInt new_size = 0;
        for (UInt i = 0; i < vect.getSize(); ++i) {
          UInt new_i = renumbering(i);
          if(new_i != UInt(-1)) {
            memcpy(tmp.storage() + new_i * nb_component,
                   vect.storage() + i *nb_component,
                   nb_component * sizeof(T));
            ++new_size;
          }
        }
        tmp.resize(new_size);
        vect.copy(tmp);
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
template<typename T, typename SupportType>
void ElementTypeMapArray<T, SupportType>::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "ElementTypeMapArray<" << debug::demangle(typeid(T).name()) << "> [" << std::endl;
  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;

    const DataMap & data = this->getData(gt);
    typename DataMap::const_iterator it;
    for(it = data.begin(); it != data.end(); ++it) {
      stream << space << space << ElementTypeMapArray::printType(it->first, gt) << " [" << std::endl;
      it->second->printself(stream, indent + 3);
      stream << space << space << " ]" << std::endl;
    }
  }
  stream << space << "]" << std::endl;
}


/* -------------------------------------------------------------------------- */
/* SupportType Iterator                                                       */
/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
ElementTypeMap<Stored, SupportType>::type_iterator::type_iterator(DataMapIterator & list_begin,
                                                                 DataMapIterator & list_end,
                                                                 UInt dim, ElementKind ek) :
  list_begin(list_begin), list_end(list_end), dim(dim), kind(ek) {
}


/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
ElementTypeMap<Stored, SupportType>::type_iterator::type_iterator(const type_iterator & it) :
  list_begin(it.list_begin), list_end(it.list_end), dim(it.dim), kind(it.kind) {
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
typename ElementTypeMap<Stored, SupportType>::type_iterator & ElementTypeMap<Stored, SupportType>::type_iterator::operator=(const type_iterator & it) {
  if(this != &it) {
    list_begin = it.list_begin;
    list_end = it.list_end;
    dim = it.dim;
    kind = it.kind;
  }
  return *this;
}
/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
inline typename ElementTypeMap<Stored, SupportType>::type_iterator::reference
ElementTypeMap<Stored, SupportType>::type_iterator::operator*() {
  return list_begin->first;
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
inline typename ElementTypeMap<Stored, SupportType>::type_iterator::reference
ElementTypeMap<Stored, SupportType>::type_iterator::operator*() const {
  return list_begin->first;
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
inline typename ElementTypeMap<Stored, SupportType>::type_iterator &
ElementTypeMap<Stored, SupportType>::type_iterator::operator++() {
  ++list_begin;
  while((list_begin != list_end) &&
	(((dim != _all_dimensions) && (dim != Mesh::getSpatialDimension(list_begin->first))) ||
	 ((kind != _ek_not_defined) && (kind != Mesh::getKind(list_begin->first)))
	 )
	)
    ++list_begin;
  return *this;
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
typename ElementTypeMap<Stored, SupportType>::type_iterator ElementTypeMap<Stored, SupportType>::type_iterator::operator++(int) {
  type_iterator tmp(*this);
  operator++();
  return tmp;
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
inline bool ElementTypeMap<Stored, SupportType>::type_iterator::operator==(const type_iterator & other) const {
  return this->list_begin == other.list_begin;
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
inline bool ElementTypeMap<Stored, SupportType>::type_iterator::operator!=(const type_iterator & other) const {
  return this->list_begin != other.list_begin;
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
inline typename ElementTypeMap<Stored, SupportType>::type_iterator
ElementTypeMap<Stored, SupportType>::firstType(UInt dim, GhostType ghost_type, ElementKind kind) const {
  typename DataMap::const_iterator b,e;
  b = getData(ghost_type).begin();
  e = getData(ghost_type).end();

  // loop until the first valid type
  while((b != e) &&
	(((dim != _all_dimensions) && (dim != Mesh::getSpatialDimension(b->first))) ||
	 ((kind != _ek_not_defined) && (kind != Mesh::getKind(b->first)))))
    ++b;

  return typename ElementTypeMap<Stored, SupportType>::type_iterator(b, e, dim, kind);
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
inline typename ElementTypeMap<Stored, SupportType>::type_iterator
ElementTypeMap<Stored, SupportType>::lastType(UInt dim, GhostType ghost_type, ElementKind kind) const {
  typename DataMap::const_iterator e;
  e = getData(ghost_type).end();
  return typename ElementTypeMap<Stored, SupportType>::type_iterator(e, e, dim, kind);
}

__END_AKANTU__

#endif /* __AKANTU_ELEMENT_TYPE_MAP_TMPL_HH__ */
