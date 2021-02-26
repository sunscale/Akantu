/**
 * @file   mesh_iterators.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Jul 16 2015
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  Set of helper classes to have fun with range based for
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
#include "aka_named_argument.hh"
#include "aka_static_if.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MESH_ITERATORS_HH_
#define AKANTU_MESH_ITERATORS_HH_

namespace akantu {

class MeshElementsByTypes {
  using elements_iterator = Array<Element>::scalar_iterator;

public:
  explicit MeshElementsByTypes(const Array<Element> & elements) {
    this->elements.copy(elements);
    std::sort(this->elements.begin(), this->elements.end());
  }

  /* ------------------------------------------------------------------------ */
  class MeshElementsRange {
  public:
    MeshElementsRange() = default;

    MeshElementsRange(const elements_iterator & begin,
                      const elements_iterator & end)
        : type((*begin).type), ghost_type((*begin).ghost_type), begin(begin),
          end(end) {}

    AKANTU_GET_MACRO(Type, type, ElementType);
    AKANTU_GET_MACRO(GhostType, ghost_type, GhostType);

    const Array<UInt> & getElements() {
      elements.resize(end - begin);
      auto el_it = elements.begin();
      for (auto it = begin; it != end; ++it, ++el_it) {
        *el_it = it->element;
      }

      return elements;
    }

  private:
    ElementType type{_not_defined};
    GhostType ghost_type{_casper};
    elements_iterator begin;
    elements_iterator end;
    Array<UInt> elements;
  };

  /* ------------------------------------------------------------------------ */
  class iterator {
    struct element_comparator {
      bool operator()(const Element & lhs, const Element & rhs) const {
        return ((rhs == ElementNull) || std::tie(lhs.ghost_type, lhs.type) <
                                            std::tie(rhs.ghost_type, rhs.type));
      }
    };

  public:
    iterator(const iterator &) = default;
    iterator(const elements_iterator & first, const elements_iterator & last)
        : range(std::equal_range(first, last, *first, element_comparator())),
          first(first), last(last) {}

    decltype(auto) operator*() const {
      return MeshElementsRange(range.first, range.second);
    }

    iterator operator++() {
      first = range.second;
      range = std::equal_range(first, last, *first, element_comparator());
      return *this;
    }

    bool operator==(const iterator & other) const {
      return (first == other.first and last == other.last);
    }

    bool operator!=(const iterator & other) const {
      return (not operator==(other));
    }

  private:
    std::pair<elements_iterator, elements_iterator> range;
    elements_iterator first;
    elements_iterator last;
  };

  iterator begin() { return iterator(elements.begin(), elements.end()); }
  iterator end() { return iterator(elements.end(), elements.end()); }

private:
  Array<Element> elements;
};

/* -------------------------------------------------------------------------- */
namespace mesh_iterators {
  namespace details {
    template <class internal_iterator> class delegated_iterator {
    public:
      using value_type = std::remove_pointer_t<
          typename internal_iterator::value_type::second_type>;
      using difference_type = std::ptrdiff_t;
      using pointer = value_type *;
      using reference = value_type &;
      using iterator_category = std::input_iterator_tag;

      explicit delegated_iterator(internal_iterator it) : it(std::move(it)) {}

      decltype(auto) operator*() {
        return std::forward<decltype(*(it->second))>(*(it->second));
      }

      delegated_iterator operator++() {
        ++it;
        return *this;
      }

      bool operator==(const delegated_iterator & other) const {
        return other.it == it;
      }

      bool operator!=(const delegated_iterator & other) const {
        return other.it != it;
      }

    private:
      internal_iterator it;
    };
  } // namespace details
} // namespace mesh_iterators

/* -------------------------------------------------------------------------- */
template <class Func>
void for_each_element(UInt nb_elements, const Array<UInt> & filter_elements,
                      Func && function) {
  if (filter_elements != empty_filter) {
    std::for_each(filter_elements.begin(), filter_elements.end(),
                  std::forward<Func>(function));
  } else {
    auto && range = arange(nb_elements);
    std::for_each(range.begin(), range.end(), std::forward<Func>(function));
  }
}

namespace {
  DECLARE_NAMED_ARGUMENT(element_filter);
}

/* -------------------------------------------------------------------------- */
template <class Func, typename... pack>
void for_each_element(const Mesh & mesh, Func && function, pack &&... _pack) {
  auto requested_ghost_type = OPTIONAL_NAMED_ARG(ghost_type, _casper);
  const ElementTypeMapArray<UInt> * filter =
      OPTIONAL_NAMED_ARG(element_filter, nullptr);

  bool all_ghost_types = requested_ghost_type == _casper;

  auto spatial_dimension =
      OPTIONAL_NAMED_ARG(spatial_dimension, mesh.getSpatialDimension());
  auto element_kind = OPTIONAL_NAMED_ARG(element_kind, _ek_not_defined);

  for (auto ghost_type : ghost_types) {
    if ((not(ghost_type == requested_ghost_type)) and (not all_ghost_types)) {
      continue;
    }

    auto element_types =
        mesh.elementTypes(spatial_dimension, ghost_type, element_kind);

    if (filter) {
      element_types =
          filter->elementTypes(spatial_dimension, ghost_type, element_kind);
    }

    for (auto type : element_types) {
      const Array<UInt> * filter_array;

      if (filter) {
        filter_array = &((*filter)(type, ghost_type));
      } else {
        filter_array = &empty_filter;
      }

      auto nb_elements = mesh.getNbElement(type, ghost_type);

      for_each_element(nb_elements, *filter_array, [&](auto && el) {
        auto element = Element{type, el, ghost_type};
        std::forward<Func>(function)(element);
      });
    }
  }
}

} // namespace akantu

#endif /* AKANTU_MESH_ITERATORS_HH_ */
