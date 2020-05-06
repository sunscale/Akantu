/**
 * @file   aka_filter_iterator.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  jeu déc 12 2019
 *
 * @brief A Documented file.
 *
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
#include <iterator>
#include <utility>
/* -------------------------------------------------------------------------- */

#ifndef AKA_FILTER_ITERATOR_H
#define AKA_FILTER_ITERATOR_H

#ifndef AKANTU_ITERATORS_NAMESPACE
#define AKANTU_ITERATORS_NAMESPACE akantu
#endif

namespace AKANTU_ITERATORS_NAMESPACE {

/* -------------------------------------------------------------------------- */
namespace iterators {
  template <class filter_iterator_t, class container_iterator_t>
  class FilterIterator {
  public:
    using value_type =
        decltype(std::declval<container_iterator_t>().operator[](0));
    using difference_type = typename filter_iterator_t::difference_type;
    using pointer = std::decay_t<value_type> *;
    using reference = value_type &;
    using iterator_category = typename filter_iterator_t::iterator_category;

    FilterIterator(filter_iterator_t filter_it,
                   container_iterator_t container_begin)
        : filter_it(std::move(filter_it)),
          container_begin(std::move(container_begin)) {}

    FilterIterator(const FilterIterator &) = default;

    FilterIterator & operator++() {
      ++filter_it;
      return *this;
    }

    decltype(auto) operator*() {
      auto container_it = this->container_begin + *this->filter_it;
      return *container_it;
    }

    decltype(auto) operator*() const {
      auto container_it = this->container_begin + *this->filter_it;
      return *container_it;
    }

    bool operator==(const FilterIterator & other) const {
      return (filter_it == other.filter_it) and
             (container_begin == other.container_begin);
    }

    bool operator!=(const FilterIterator & other) const {
      return filter_it != other.filter_it;
    }

  private:
    filter_iterator_t filter_it;
    container_iterator_t container_begin;
  };

  template <class filter_iterator_t, class container_iterator_t>
  decltype(auto) make_filter_iterator(filter_iterator_t && filter_it,
                                      container_iterator_t && container_begin) {
    return FilterIterator<filter_iterator_t, container_iterator_t>(
        std::forward<filter_iterator_t>(filter_it),
        std::forward<container_iterator_t>(container_begin));
  }
} // namespace iterators

namespace containers {
  template <class filter_t, class Container> class FilterAdaptor {
  public:
    FilterAdaptor(filter_t && filter, Container && container)
        : filter(std::forward<filter_t>(filter)),
          container(std::forward<Container>(container)) {
      static_assert(
          std::is_same<typename decltype(container.begin())::iterator_category,
                       std::random_access_iterator_tag>::value,
          "Containers must all have random iterators");
    }

    decltype(auto) begin() const {
      return iterators::make_filter_iterator(filter.begin(), container.begin());
    }
    decltype(auto) begin() {
      return iterators::make_filter_iterator(filter.begin(), container.begin());
    }

    decltype(auto) end() const {
      return iterators::make_filter_iterator(filter.end(), container.begin());
    }
    decltype(auto) end() {
      return iterators::make_filter_iterator(filter.end(), container.begin());
    }

  private:
    filter_t filter;
    Container container;
  };
} // namespace containers

template <class filter_t, class Container>
decltype(auto) filter(filter_t && filter, Container && container) {
  return containers::FilterAdaptor<filter_t, Container>(
      std::forward<filter_t>(filter), std::forward<Container>(container));
}

} // namespace AKANTU_ITERATORS_NAMESPACE

namespace std {
template <class filter_iterator_t, class container_iterator_t>
struct iterator_traits<::AKANTU_ITERATORS_NAMESPACE::iterators::FilterIterator<
    filter_iterator_t, container_iterator_t>> {
private:
  using iterator_type =
      typename ::AKANTU_ITERATORS_NAMESPACE::iterators::FilterIterator<
          filter_iterator_t, container_iterator_t>;

public:
  using iterator_category = typename iterator_type::iterator_category;
  using value_type = typename iterator_type::value_type;
  using difference_type = typename iterator_type::difference_type;
  using pointer = typename iterator_type::pointer;
  using reference = typename iterator_type::reference;
};
} // namespace std

#endif // AKA_FILTER_ITERATOR_H
