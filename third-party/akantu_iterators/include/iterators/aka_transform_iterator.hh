/**
 * @file   aka_transform_iterator.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  jeu déc 12 2019
 *
 * @brief transform adaptors
 *
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * akantu-iterators is free  software: you can redistribute it and/or  modify it
 * under the terms  of the  GNU Lesser  General Public  License as  published by
 * the Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * akantu-iterators is  distributed in the  hope that it  will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public
 * License  for more details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with akantu-iterators. If not, see <http://www.gnu.org/licenses/>.
 *
 */
/* -------------------------------------------------------------------------- */
#include "aka_compatibilty_with_cpp_standard.hh"
/* -------------------------------------------------------------------------- */
#include <utility>
/* -------------------------------------------------------------------------- */

#ifndef AKA_TRANSFORM_ITERATOR_H
#define AKA_TRANSFORM_ITERATOR_H

#ifndef AKANTU_ITERATORS_NAMESPACE
#define AKANTU_ITERATORS_NAMESPACE akantu
#endif

namespace AKANTU_ITERATORS_NAMESPACE {

namespace iterators {
  template <class iterator_t, class operator_t>
  class transform_adaptor_iterator {
  public:
    using value_type = decltype(std::declval<operator_t>()(
        std::declval<typename iterator_t::value_type>()));
    using difference_type = typename iterator_t::difference_type;
    using pointer = std::decay_t<value_type> *;
    using reference = value_type &;
    using iterator_category = typename iterator_t::iterator_category;

    transform_adaptor_iterator(iterator_t it, operator_t op)
        : it(std::move(it)), op(op) {}
    transform_adaptor_iterator(const transform_adaptor_iterator &) = default;

    transform_adaptor_iterator & operator++() {
      ++it;
      return *this;
    }

    decltype(auto) operator*() { return op(std::forward<decltype(*it)>(*it)); }

    bool operator==(const transform_adaptor_iterator & other) const {
      return (it == other.it);
    }

    bool operator!=(const transform_adaptor_iterator & other) const {
      return not operator==(other);
    }

    template <class iterator_category_ = iterator_category,
              std::enable_if_t<aka::is_iterator_category_at_least<
                  iterator_category_,
                  std::random_access_iterator_tag>::value> * = nullptr>
    difference_type operator-(const transform_adaptor_iterator & other) {
      return other - *this;
    }

  private:
    iterator_t it;
    operator_t op;
  };

  template <class iterator_t, class operator_t>
  decltype(auto) make_transform_adaptor_iterator(iterator_t it, operator_t op) {
    return transform_adaptor_iterator<iterator_t, operator_t>(
        it, std::forward<operator_t>(op));
  }

} // namespace iterators

namespace containers {
  template <class container_t, class operator_t>
  class TransformIteratorAdaptor {
  public:
    // using const_iterator = typename
    // std::decay_t<container_t>::const_iterator; using iterator = typename
    // std::decay_t<container_t>::iterator;

    TransformIteratorAdaptor(container_t && cont, operator_t op)
        : cont(std::forward<container_t>(cont)),
          op(std::forward<operator_t>(op)) {}

    decltype(auto) begin() const {
      return iterators::make_transform_adaptor_iterator(cont.begin(), op);
    }
    decltype(auto) begin() {
      return iterators::make_transform_adaptor_iterator(cont.begin(), op);
    }

    decltype(auto) end() const {
      return iterators::make_transform_adaptor_iterator(cont.end(), op);
    }
    decltype(auto) end() {
      return iterators::make_transform_adaptor_iterator(cont.end(), op);
    }

  private:
    container_t cont;
    operator_t op;
  };
} // namespace containers

template <class container_t, class operator_t>
decltype(auto) make_transform_adaptor(container_t && cont, operator_t && op) {
  return containers::TransformIteratorAdaptor<container_t, operator_t>(
      std::forward<container_t>(cont), std::forward<operator_t>(op));
}

template <class container_t>
decltype(auto) make_keys_adaptor(container_t && cont) {
  return make_transform_adaptor(
      std::forward<container_t>(cont),
      [](auto && pair) -> const auto & { return pair.first; });
}

template <class container_t>
decltype(auto) make_values_adaptor(container_t && cont) {
  return make_transform_adaptor(
      std::forward<container_t>(cont),
      [](auto && pair) -> auto & { return pair.second; });
}

template <class container_t>
decltype(auto) make_dereference_adaptor(container_t && cont) {
  return make_transform_adaptor(
      std::forward<container_t>(cont),
      [](auto && value) -> decltype(*value) { return *value; });
}

} // namespace AKANTU_ITERATORS_NAMESPACE

#endif // AKA_TRANSFORM_ITERATOR_H
