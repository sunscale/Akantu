/**
 * @file   element_type_map_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Aug 31 2011
 * @date last modification: Fri Oct 02 2015
 *
 * @brief  implementation of template functions of the ElementTypeMap and
 * ElementTypeMapArray classes
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "aka_static_if.hh"
#include "element_type_map.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */
#include "element_type_conversion.hh"
/* -------------------------------------------------------------------------- */
#include <tuple>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_ELEMENT_TYPE_MAP_TMPL_HH__
#define __AKANTU_ELEMENT_TYPE_MAP_TMPL_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
/* ElementTypeMap                                                             */
/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
inline std::string
ElementTypeMap<Stored, SupportType>::printType(const SupportType & type,
                                               const GhostType & ghost_type) {
  std::stringstream sstr;
  sstr << "(" << ghost_type << ":" << type << ")";
  return sstr.str();
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
inline bool ElementTypeMap<Stored, SupportType>::exists(
    const SupportType & type, const GhostType & ghost_type) const {
  return this->getData(ghost_type).find(type) !=
         this->getData(ghost_type).end();
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
inline const Stored & ElementTypeMap<Stored, SupportType>::
operator()(const SupportType & type, const GhostType & ghost_type) const {
  auto it = this->getData(ghost_type).find(type);

  if (it == this->getData(ghost_type).end())
    AKANTU_SILENT_EXCEPTION("No element of type "
                            << ElementTypeMap::printType(type, ghost_type)
                            << " in this ElementTypeMap<"
                            << debug::demangle(typeid(Stored).name())
                            << "> class");
  return it->second;
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
inline Stored & ElementTypeMap<Stored, SupportType>::
operator()(const SupportType & type, const GhostType & ghost_type) {
  return this->getData(ghost_type)[type];
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
inline Stored & ElementTypeMap<Stored, SupportType>::
operator()(const Stored & insert, const SupportType & type,
           const GhostType & ghost_type) {
  auto it = this->getData(ghost_type).find(type);

  if (it != this->getData(ghost_type).end()) {
    AKANTU_SILENT_EXCEPTION("Element of type "
                            << ElementTypeMap::printType(type, ghost_type)
                            << " already in this ElementTypeMap<"
                            << debug::demangle(typeid(Stored).name())
                            << "> class");
  } else {
    auto & data = this->getData(ghost_type);
    const auto & res =
        data.insert(std::pair<ElementType, Stored>(type, insert));
    it = res.first;
  }

  return it->second;
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
inline typename ElementTypeMap<Stored, SupportType>::DataMap &
ElementTypeMap<Stored, SupportType>::getData(GhostType ghost_type) {
  if (ghost_type == _not_ghost)
    return data;

  return ghost_data;
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
inline const typename ElementTypeMap<Stored, SupportType>::DataMap &
ElementTypeMap<Stored, SupportType>::getData(GhostType ghost_type) const {
  if (ghost_type == _not_ghost)
    return data;

  return ghost_data;
}

/* -------------------------------------------------------------------------- */
/// Works only if stored is a pointer to a class with a printself method
template <class Stored, typename SupportType>
void ElementTypeMap<Stored, SupportType>::printself(std::ostream & stream,
                                                    int indent) const {
  std::string space;
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
    ;

  stream << space << "ElementTypeMap<" << debug::demangle(typeid(Stored).name())
         << "> [" << std::endl;
  for (auto gt : ghost_types) {
    const DataMap & data = getData(gt);
    for (auto & pair : data) {
      stream << space << space << ElementTypeMap::printType(pair.first, gt)
             << std::endl;
    }
  }

  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
ElementTypeMap<Stored, SupportType>::ElementTypeMap() = default;

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
ElementTypeMap<Stored, SupportType>::~ElementTypeMap() = default;

/* -------------------------------------------------------------------------- */
/* ElementTypeMapArray                                                        */
/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType>
inline Array<T> & ElementTypeMapArray<T, SupportType>::alloc(
    UInt size, UInt nb_component, const SupportType & type,
    const GhostType & ghost_type, const T & default_value) {
  std::string ghost_id = "";
  if (ghost_type == _ghost)
    ghost_id = ":ghost";

  Array<T> * tmp;

  auto it = this->getData(ghost_type).find(type);

  if (it == this->getData(ghost_type).end()) {
    std::stringstream sstr;
    sstr << this->id << ":" << type << ghost_id;
    tmp = &(Memory::alloc<T>(sstr.str(), size, nb_component, default_value));
    std::stringstream sstrg;
    sstrg << ghost_type;
    // tmp->setTag(sstrg.str());
    this->getData(ghost_type)[type] = tmp;
  } else {
    AKANTU_DEBUG_INFO(
        "The vector "
        << this->id << this->printType(type, ghost_type)
        << " already exists, it is resized instead of allocated.");
    tmp = it->second;
    it->second->resize(size);
  }

  return *tmp;
}

/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType>
inline void
ElementTypeMapArray<T, SupportType>::alloc(UInt size, UInt nb_component,
                                           const SupportType & type,
                                           const T & default_value) {
  this->alloc(size, nb_component, type, _not_ghost, default_value);
  this->alloc(size, nb_component, type, _ghost, default_value);
}

/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType>
inline void ElementTypeMapArray<T, SupportType>::free() {
  AKANTU_DEBUG_IN();

  for (auto && gt : ghost_types) {
    auto & data = this->getData(gt);
    for (auto & pair : data) {
      dealloc(pair.second->getID());
    }
    data.clear();
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType>
inline void ElementTypeMapArray<T, SupportType>::clear() {
  for (auto gt : ghost_types) {
    auto & data = this->getData(gt);
    for (auto & vect : data) {
      vect.second->clear();
    }
  }
}

/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType>
inline const Array<T> & ElementTypeMapArray<T, SupportType>::
operator()(const SupportType & type, const GhostType & ghost_type) const {
  auto it = this->getData(ghost_type).find(type);

  if (it == this->getData(ghost_type).end())
    AKANTU_SILENT_EXCEPTION("No element of type "
                            << ElementTypeMapArray::printType(type, ghost_type)
                            << " in this const ElementTypeMapArray<"
                            << debug::demangle(typeid(T).name()) << "> class(\""
                            << this->id << "\")");
  return *(it->second);
}

/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType>
inline Array<T> & ElementTypeMapArray<T, SupportType>::
operator()(const SupportType & type, const GhostType & ghost_type) {
  auto it = this->getData(ghost_type).find(type);

  if (it == this->getData(ghost_type).end())
    AKANTU_SILENT_EXCEPTION("No element of type "
                            << ElementTypeMapArray::printType(type, ghost_type)
                            << " in this ElementTypeMapArray<"
                            << debug::demangle(typeid(T).name())
                            << "> class (\"" << this->id << "\")");

  return *(it->second);
}

/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType>
inline void
ElementTypeMapArray<T, SupportType>::setArray(const SupportType & type,
                                              const GhostType & ghost_type,
                                              const Array<T> & vect) {
  auto it = this->getData(ghost_type).find(type);

  if (AKANTU_DEBUG_TEST(dblWarning) && it != this->getData(ghost_type).end() &&
      it->second != &vect) {
    AKANTU_DEBUG_WARNING(
        "The Array "
        << this->printType(type, ghost_type)
        << " is already registred, this call can lead to a memory leak.");
  }

  this->getData(ghost_type)[type] = &(const_cast<Array<T> &>(vect));
}

/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType>
inline void ElementTypeMapArray<T, SupportType>::onElementsRemoved(
    const ElementTypeMapArray<UInt> & new_numbering) {
  for (auto gt : ghost_types) {
    for (auto & type :
         new_numbering.elementTypes(_all_dimensions, gt, _ek_not_defined)) {
      SupportType support_type = convertType<ElementType, SupportType>(type);
      if (this->exists(support_type, gt)) {
        const auto & renumbering = new_numbering(type, gt);
        if (renumbering.size() == 0)
          continue;

        auto & vect = this->operator()(support_type, gt);
        auto nb_component = vect.getNbComponent();
        Array<T> tmp(renumbering.size(), nb_component);
        UInt new_size = 0;

        for (UInt i = 0; i < vect.size(); ++i) {
          UInt new_i = renumbering(i);
          if (new_i != UInt(-1)) {
            memcpy(tmp.storage() + new_i * nb_component,
                   vect.storage() + i * nb_component, nb_component * sizeof(T));
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
template <typename T, typename SupportType>
void ElementTypeMapArray<T, SupportType>::printself(std::ostream & stream,
                                                    int indent) const {
  std::string space;
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
    ;

  stream << space << "ElementTypeMapArray<" << debug::demangle(typeid(T).name())
         << "> [" << std::endl;
  for (UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType)g;

    const DataMap & data = this->getData(gt);
    typename DataMap::const_iterator it;
    for (it = data.begin(); it != data.end(); ++it) {
      stream << space << space << ElementTypeMapArray::printType(it->first, gt)
             << " [" << std::endl;
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
ElementTypeMap<Stored, SupportType>::type_iterator::type_iterator(
    DataMapIterator & list_begin, DataMapIterator & list_end, UInt dim,
    ElementKind ek)
    : list_begin(list_begin), list_end(list_end), dim(dim), kind(ek) {}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
ElementTypeMap<Stored, SupportType>::type_iterator::type_iterator(
    const type_iterator & it)
    : list_begin(it.list_begin), list_end(it.list_end), dim(it.dim),
      kind(it.kind) {}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
typename ElementTypeMap<Stored, SupportType>::type_iterator &
ElementTypeMap<Stored, SupportType>::type_iterator::
operator=(const type_iterator & it) {
  if (this != &it) {
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
  while ((list_begin != list_end) &&
         (((dim != _all_dimensions) &&
           (dim != Mesh::getSpatialDimension(list_begin->first))) ||
          ((kind != _ek_not_defined) &&
           (kind != Mesh::getKind(list_begin->first)))))
    ++list_begin;
  return *this;
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
typename ElementTypeMap<Stored, SupportType>::type_iterator
    ElementTypeMap<Stored, SupportType>::type_iterator::operator++(int) {
  type_iterator tmp(*this);
  operator++();
  return tmp;
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
inline bool ElementTypeMap<Stored, SupportType>::type_iterator::
operator==(const type_iterator & other) const {
  return this->list_begin == other.list_begin;
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
inline bool ElementTypeMap<Stored, SupportType>::type_iterator::
operator!=(const type_iterator & other) const {
  return this->list_begin != other.list_begin;
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
typename ElementTypeMap<Stored, SupportType>::ElementTypesIteratorHelper
ElementTypeMap<Stored, SupportType>::elementTypesImpl(
    UInt dim, GhostType ghost_type, ElementKind kind) const {
  return ElementTypesIteratorHelper(*this, dim, ghost_type, kind);
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
template <typename... pack>
typename ElementTypeMap<Stored, SupportType>::ElementTypesIteratorHelper
ElementTypeMap<Stored, SupportType>::elementTypesImpl(
    const use_named_args_t & /*unused*/, pack &&... _pack) const {
  return ElementTypesIteratorHelper(*this, use_named_args, _pack...);
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
template <typename... pack>
decltype(auto)
ElementTypeMap<Stored, SupportType>::elementTypes(pack &&... _pack) const {
  auto && first_arg = std::get<0>(std::forward_as_tuple(_pack...));

  return static_if(is_named_argument<std::decay_t<decltype(first_arg)>>{})
      .then([&](auto && /*a*/) {
        return elementTypesImpl(use_named_args,
                                    std::forward<decltype(_pack)>(_pack)...);
      })
      .else_([&](auto && /*a*/) {
        return elementTypesImpl(std::forward<decltype(_pack)>(_pack)...);
      })(std::forward<decltype(first_arg)>(first_arg));
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
inline typename ElementTypeMap<Stored, SupportType>::type_iterator
ElementTypeMap<Stored, SupportType>::firstType(UInt dim, GhostType ghost_type,
                                               ElementKind kind) const {
  typename DataMap::const_iterator b, e;
  b = getData(ghost_type).begin();
  e = getData(ghost_type).end();

  // loop until the first valid type
  while ((b != e) &&
         (((dim != _all_dimensions) &&
           (dim != Mesh::getSpatialDimension(b->first))) ||
          ((kind != _ek_not_defined) && (kind != Mesh::getKind(b->first)))))
    ++b;

  return typename ElementTypeMap<Stored, SupportType>::type_iterator(b, e, dim,
                                                                     kind);
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
inline typename ElementTypeMap<Stored, SupportType>::type_iterator
ElementTypeMap<Stored, SupportType>::lastType(UInt dim, GhostType ghost_type,
                                              ElementKind kind) const {
  typename DataMap::const_iterator e;
  e = getData(ghost_type).end();
  return typename ElementTypeMap<Stored, SupportType>::type_iterator(e, e, dim,
                                                                     kind);
}

/* -------------------------------------------------------------------------- */

/// standard output stream operator
template <class Stored, typename SupportType>
inline std::ostream &
operator<<(std::ostream & stream,
           const ElementTypeMap<Stored, SupportType> & _this) {
  _this.printself(stream);
  return stream;
}

/* -------------------------------------------------------------------------- */
class ElementTypeMapArrayInializer {
public:
  ElementTypeMapArrayInializer(UInt spatial_dimension = _all_dimensions,
                               UInt nb_component = 1,
                               const GhostType & ghost_type = _not_ghost,
                               const ElementKind & element_kind = _ek_regular)
      : spatial_dimension(spatial_dimension), nb_component(nb_component),
        ghost_type(ghost_type), element_kind(element_kind) {}

  const GhostType & ghostType() const { return ghost_type; }

protected:
  UInt spatial_dimension;
  UInt nb_component;
  GhostType ghost_type;
  ElementKind element_kind;
};

/* -------------------------------------------------------------------------- */
class MeshElementTypeMapArrayInializer : public ElementTypeMapArrayInializer {
public:
  MeshElementTypeMapArrayInializer(
      const Mesh & mesh, UInt nb_component = 1,
      UInt spatial_dimension = _all_dimensions,
      const GhostType & ghost_type = _not_ghost,
      const ElementKind & element_kind = _ek_regular,
      bool with_nb_element = false, bool with_nb_nodes_per_element = false)
      : ElementTypeMapArrayInializer(spatial_dimension, nb_component,
                                     ghost_type, element_kind),
        mesh(mesh), with_nb_element(with_nb_element),
        with_nb_nodes_per_element(with_nb_nodes_per_element) {}

  decltype(auto) elementTypes() const {
    return mesh.elementTypes(this->spatial_dimension, this->ghost_type,
                             this->element_kind);
  }

  virtual UInt size(const ElementType & type) const {
    if (with_nb_element)
      return mesh.getNbElement(type, this->ghost_type);

    return 0;
  }

  UInt nbComponent(const ElementType & type) const {
    if (with_nb_nodes_per_element)
      return (this->nb_component * mesh.getNbNodesPerElement(type));

    return this->nb_component;
  }

protected:
  const Mesh & mesh;
  bool with_nb_element;
  bool with_nb_nodes_per_element;
};

/* -------------------------------------------------------------------------- */
class FEEngineElementTypeMapArrayInializer
    : public MeshElementTypeMapArrayInializer {
public:
  FEEngineElementTypeMapArrayInializer(
      const FEEngine & fe_engine, UInt nb_component = 1,
      UInt spatial_dimension = _all_dimensions,
      const GhostType & ghost_type = _not_ghost,
      const ElementKind & element_kind = _ek_regular);

  UInt size(const ElementType & type) const override;

  using ElementTypesIteratorHelper =
      ElementTypeMapArray<Real, ElementType>::ElementTypesIteratorHelper;

  ElementTypesIteratorHelper elementTypes() const;

protected:
  const FEEngine & fe_engine;
};

/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType>
template <class Func>
void ElementTypeMapArray<T, SupportType>::initialize(const Func & f,
                                                     const T & default_value) {
  for (auto & type : f.elementTypes()) {
    if (not this->exists(type, f.ghostType()))
      this->alloc(f.size(type), f.nbComponent(type), type, f.ghostType(),
                  default_value);
  }
}

/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType>
template <typename... pack>
void ElementTypeMapArray<T, SupportType>::initialize(const Mesh & mesh,
                                                     pack &&... _pack) {
  GhostType requested_ghost_type = OPTIONAL_NAMED_ARG(ghost_type, _casper);
  bool all_ghost_types = requested_ghost_type == _casper;

  for (auto ghost_type : ghost_types) {
    if ((not(ghost_type == requested_ghost_type)) and (not all_ghost_types))
      continue;

    this->initialize(
        MeshElementTypeMapArrayInializer(
            mesh, OPTIONAL_NAMED_ARG(nb_component, 1),
            OPTIONAL_NAMED_ARG(spatial_dimension, mesh.getSpatialDimension()),
            ghost_type, OPTIONAL_NAMED_ARG(element_kind, _ek_regular),
            OPTIONAL_NAMED_ARG(with_nb_element, false),
            OPTIONAL_NAMED_ARG(with_nb_nodes_per_element, false)),
        OPTIONAL_NAMED_ARG(default_value, T()));
  }
}

/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType>
template <typename... pack>
void ElementTypeMapArray<T, SupportType>::initialize(const FEEngine & fe_engine,
                                                     pack &&... _pack) {
  bool all_ghost_types = OPTIONAL_NAMED_ARG(all_ghost_types, true);
  GhostType requested_ghost_type = OPTIONAL_NAMED_ARG(ghost_type, _not_ghost);

  for (auto ghost_type : ghost_types) {
    if ((not(ghost_type == requested_ghost_type)) and (not all_ghost_types))
      continue;

    this->initialize(FEEngineElementTypeMapArrayInializer(
                         fe_engine, OPTIONAL_NAMED_ARG(nb_component, 1),
                         OPTIONAL_NAMED_ARG(spatial_dimension, UInt(-2)),
                         ghost_type,
                         OPTIONAL_NAMED_ARG(element_kind, _ek_regular)),
                     OPTIONAL_NAMED_ARG(default_value, T()));
  }
}

/* -------------------------------------------------------------------------- */
template <class T, typename SupportType>
inline T & ElementTypeMapArray<T, SupportType>::
operator()(const Element & element) {
  return this->operator()(element.type, element.ghost_type)(element.element);
}

/* -------------------------------------------------------------------------- */
template <class T, typename SupportType>
inline const T & ElementTypeMapArray<T, SupportType>::
operator()(const Element & element) const {
  return this->operator()(element.type, element.ghost_type)(element.element);
}

} // namespace akantu

#endif /* __AKANTU_ELEMENT_TYPE_MAP_TMPL_HH__ */
