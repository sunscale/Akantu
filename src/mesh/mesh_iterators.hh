/**
 * @file   mesh_iterators.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  Wed Aug 09 2017
 *
 * @brief Set of helper classes to have fun with range based for
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "mesh.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MESH_ITERATORS_HH__
#define __AKANTU_MESH_ITERATORS_HH__

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

    AKANTU_GET_MACRO(Type, type, const ElementType &);
    AKANTU_GET_MACRO(GhostType, ghost_type, const GhostType &);

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
        bool res =
            ((lhs.kind < rhs.kind) ||
             ((lhs.kind == rhs.kind) && ((lhs.ghost_type < rhs.ghost_type) ||
                                         ((lhs.ghost_type == rhs.ghost_type) &&
                                          ((lhs.type < rhs.type))))));
        return res;
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
      explicit delegated_iterator(internal_iterator it) : it(std::move(it)) {}

      decltype(auto) operator*() {
        return std::forward<decltype(*(it->second))>(*(it->second));
      }

      delegated_iterator operator++() {
        ++it;
        return *this;
      }

      bool operator==(const delegated_iterator & other) {
        return other.it == it;
      }

      bool operator!=(const delegated_iterator & other) {
        return other.it != it;
      }

    private:
      internal_iterator it;
    };

    template <class M> class MeshElementGroups {
    public:
      explicit MeshElementGroups(M && mesh) : mesh(std::forward<M>(mesh)) {}

      decltype(auto) begin() {
        return delegated_iterator<decltype(mesh.element_group_begin())>(
            mesh.element_group_begin());
      }

      decltype(auto) begin() const {
        return delegated_iterator<decltype(mesh.element_group_begin())>(
            mesh.element_group_begin());
      }

      decltype(auto) end() {
        return delegated_iterator<decltype(mesh.element_group_end())>(
            mesh.element_group_end());
      }

      decltype(auto) end() const {
        return delegated_iterator<decltype(mesh.element_group_end())>(
            mesh.element_group_end());
      }

    private:
      M && mesh;
    };

    template <class M> class MeshNodeGroups {
    public:
      explicit MeshNodeGroups(M && mesh) : mesh(std::forward<M>(mesh)) {}

      decltype(auto) begin() {
        return delegated_iterator<decltype(mesh.node_group_begin())>(
            mesh.node_group_begin());
      }

      decltype(auto) begin() const {
        return delegated_iterator<decltype(mesh.node_group_begin())>(
            mesh.node_group_begin());
      }

      decltype(auto) end() {
        return delegated_iterator<decltype(mesh.node_group_end())>(
            mesh.node_group_end());
      }

      decltype(auto) end() const {
        return delegated_iterator<decltype(mesh.node_group_end())>(
            mesh.node_group_end());
      }

    private:
      M && mesh;
    };
  } // namespace details
} // namespace mesh_iterators

template <class Mesh> decltype(auto) MeshElementGroups(Mesh && mesh) {
  return mesh_iterators::details::MeshElementGroups<Mesh>(mesh);
}

template <class Mesh> decltype(auto) MeshNodeGroups(Mesh && mesh) {
  return mesh_iterators::details::MeshNodeGroups<Mesh>(mesh);
}

} // namespace akantu

#endif /* __AKANTU_MESH_ITERATORS_HH__ */
