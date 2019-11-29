/**
 * @file   aka_iterators.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Aug 11 2017
 * @date last modification: Mon Jan 29 2018
 *
 * @brief  iterator interfaces
 *
 * @section LICENSE
 *
 * Copyright (©) 2016-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_compatibilty_with_cpp_standard.hh"
#include "aka_tuple_tools.hh"
/* -------------------------------------------------------------------------- */
#include <iterator>
#include <tuple>
#include <utility>

/* -------------------------------------------------------------------------- */
#ifndef AKANTU_AKA_ITERATORS_HH
#define AKANTU_AKA_ITERATORS_HH

#ifndef AKANTU_ITERATOR_NAMESPACE
#define AKANTU_ITERATOR_NAMESPACE akantu
#endif

namespace AKANTU_ITERATOR_NAMESPACE {

/* -------------------------------------------------------------------------- */
namespace iterators {
namespace details {
template <typename cat1, typename cat2>
using is_iterator_category_at_least =
    std::is_same<std::common_type_t<cat1, cat2>, cat2>;
}

template <class... Iterators> class ZipIterator {
private:
  using tuple_t = std::tuple<Iterators...>;

public:
  using value_type =
      std::tuple<typename std::iterator_traits<Iterators>::value_type...>;
  using difference_type = std::common_type_t<
    typename std::iterator_traits<Iterators>::difference_type...>;
  using pointer =
      std::tuple<typename std::iterator_traits<Iterators>::pointer...>;
  using reference =
      std::tuple<typename std::iterator_traits<Iterators>::reference...>;
  using iterator_category = // std::input_iterator_tag;
      std::common_type_t<
    typename std::iterator_traits<Iterators>::iterator_category...>;

  using nb_iterators = std::tuple_size<tuple_t>;

public:
  explicit ZipIterator(tuple_t iterators) : iterators(std::move(iterators)) {}

  ZipIterator(const ZipIterator & other) = default;
  ZipIterator(ZipIterator && other) noexcept = default;

  ZipIterator & operator=(const ZipIterator & other) = default;
  ZipIterator & operator=(ZipIterator && other) noexcept = default;

  /* ---------------------------------------------------------------------- */
  template <class iterator_category_ = iterator_category,
            std::enable_if_t<details::is_iterator_category_at_least<
                               iterator_category_,
                               std::bidirectional_iterator_tag>::value> * = nullptr>
  ZipIterator & operator--() {
    tuple::foreach ([](auto && it) { --it; }, iterators);
    return *this;
  }

  template <class iterator_category_ = iterator_category,
            std::enable_if_t<details::is_iterator_category_at_least<
                               iterator_category_,
                               std::bidirectional_iterator_tag>::value> * = nullptr>
  ZipIterator operator--(int a) {
    auto cpy = *this;
    this->operator--();
    return cpy;
  }

  // input iterator ++it
  ZipIterator & operator++() {
    tuple::foreach ([](auto && it) { ++it; }, iterators);
    return *this;
  }

  // input iterator it++
  ZipIterator operator++(int) {
    auto cpy = *this;
    this->operator++();
    return cpy;
  }

  // input iterator it != other_it
  bool operator!=(const ZipIterator & other) const {
    // return tuple::are_not_equal(iterators, other.iterators);
    return std::get<0>(iterators) !=
        std::get<0>(other.iterators); // helps the compiler to optimize
  }

  // input iterator dereference *it
  decltype(auto) operator*() {
    return tuple::transform([](auto && it) -> decltype(auto) { return *it; },
                            iterators);
  }

  template <class iterator_category_ = iterator_category,
            std::enable_if_t<details::is_iterator_category_at_least<
                               iterator_category_,
                               std::random_access_iterator_tag>::value> * = nullptr>
  difference_type operator-(const ZipIterator & other) {
    return std::get<0>(this->iterators) - std::get<0>(other.iterators);
  }

  // random iterator it[idx]
  template <class iterator_category_ = iterator_category,
            std::enable_if_t<details::is_iterator_category_at_least<
                               iterator_category_,
                               std::random_access_iterator_tag>::value> * = nullptr>
  decltype(auto) operator[](std::size_t idx) {
    return tuple::transform(
        [idx](auto && it) -> decltype(auto) { return it[idx]; }, iterators);
  }

  // random iterator it + n
  template <class iterator_category_ = iterator_category,
            std::enable_if_t<details::is_iterator_category_at_least<
                               iterator_category_,
                               std::random_access_iterator_tag>::value> * = nullptr>
  decltype(auto) operator+(std::size_t n) {
    return ZipIterator(std::forward<tuple_t>(tuple::transform(
                                                 [n](auto && it) -> decltype(auto) { return it + n; }, iterators)));
  }

  // random iterator it - n
  template <class iterator_category_ = iterator_category,
            std::enable_if_t<details::is_iterator_category_at_least<
                               iterator_category_,
                               std::random_access_iterator_tag>::value> * = nullptr>
  decltype(auto) operator-(std::size_t n) {
    return ZipIterator(std::forward<tuple_t>(tuple::transform(
                                                 [n](auto && it) -> decltype(auto) { return it - n; }, iterators)));
  }

  template <
    class iterator_category_ = iterator_category,
    std::enable_if_t<details::is_iterator_category_at_least<
                       iterator_category_, std::forward_iterator_tag>::value> * = nullptr>
  bool operator==(const ZipIterator & other) const {
    return not tuple::are_not_equal(iterators, other.iterators);
  }

private:
  tuple_t iterators;
};
} // namespace iterators

/* -------------------------------------------------------------------------- */
template <class... Iterators>
decltype(auto) zip_iterator(std::tuple<Iterators...> && iterators_tuple) {
  auto zip = iterators::ZipIterator<Iterators...>(
      std::forward<decltype(iterators_tuple)>(iterators_tuple));
  return zip;
}

/* -------------------------------------------------------------------------- */
namespace containers {
template <class... Containers> class ZipContainer {
  using containers_t = std::tuple<Containers...>;

public:
  explicit ZipContainer(Containers &&... containers)
      : containers(std::forward<Containers>(containers)...) {}

  decltype(auto) begin() const {
    return zip_iterator(
        tuple::transform([](auto && c) { return c.begin(); },
                         std::forward<containers_t>(containers)));
  }

  decltype(auto) end() const {
    return zip_iterator(
        tuple::transform([](auto && c) { return c.end(); },
                         std::forward<containers_t>(containers)));
  }

  decltype(auto) begin() {
    return zip_iterator(
        tuple::transform([](auto && c) { return c.begin(); },
                         std::forward<containers_t>(containers)));
  }

  decltype(auto) end() {
    return zip_iterator(
        tuple::transform([](auto && c) { return c.end(); },
                         std::forward<containers_t>(containers)));
  }

  // template <class Container = std::tuple_element<0, containers_t>,
  //           std::enable_if_t<std::is_integral<decltype(
  //               std::declval<Container>().size())>::value> * = nullptr>
  // decltype(auto) size() {
  //   return std::forward<Container>(std::get<0>(containers)).size();
  // }

private:
  containers_t containers;
};

template <class Iterator> class Range {
public:
  using iterator = Iterator;
  // ugly trick
  using const_iterator = Iterator;

  explicit Range(Iterator && it1, Iterator && it2)
      : iterators(std::forward<Iterator>(it1), std::forward<Iterator>(it2)) {}

  decltype(auto) begin() const { return std::get<0>(iterators); }
  decltype(auto) begin() { return std::get<0>(iterators); }

  decltype(auto) end() const { return std::get<1>(iterators); }
  decltype(auto) end() { return std::get<1>(iterators); }

private:
  std::tuple<Iterator, Iterator> iterators;
};
} // namespace containers

/* -------------------------------------------------------------------------- */
template <class... Containers> decltype(auto) zip(Containers &&... conts) {
  return containers::ZipContainer<Containers...>(
      std::forward<Containers>(conts)...);
}

template <class Iterator>
decltype(auto) range(Iterator && it1, Iterator && it2) {
  return containers::Range<Iterator>(std::forward<Iterator>(it1),
                                     std::forward<Iterator>(it2));
}
/* -------------------------------------------------------------------------- */
/* Arange                                                                     */
/* -------------------------------------------------------------------------- */
namespace iterators {
template <class T> class ArangeIterator {
public:
  using value_type = T;
  using pointer = T *;
  using reference = T &;
  using difference_type = size_t;
  using iterator_category = std::forward_iterator_tag;

  constexpr ArangeIterator(T value, T step) : value(value), step(step) {}
  constexpr ArangeIterator(const ArangeIterator &) = default;

  constexpr ArangeIterator & operator++() {
    value += step;
    return *this;
  }

  constexpr T operator*() const { return value; }

  constexpr bool operator==(const ArangeIterator & other) const {
    return (value == other.value) and (step == other.step);
  }

  constexpr bool operator!=(const ArangeIterator & other) const {
    return not operator==(other);
  }

private:
  T value{0};
  const T step{1};
};
} // namespace iterators

namespace containers {
template <class T> class ArangeContainer {
public:
  using iterator = iterators::ArangeIterator<T>;
  using const_iterator = iterators::ArangeIterator<T>;

  constexpr ArangeContainer(T start, T stop, T step = 1)
      : start(start), stop((stop - start) % step == 0
                           ? stop
                           : start + (1 + (stop - start) / step) * step),
        step(step) {}
  explicit constexpr ArangeContainer(T stop) : ArangeContainer(0, stop, 1) {}

  constexpr T operator[](size_t i) {
    T val = start + i * step;
    assert(val < stop && "i is out of range");
    return val;
  }

  constexpr T size() { return (stop - start) / step; }

  constexpr iterator begin() { return iterator(start, step); }
  constexpr iterator end() { return iterator(stop, step); }

private:
  const T start{0}, stop{0}, step{1};
};
} // namespace containers

template <class T,
          typename = std::enable_if_t<std::is_integral<std::decay_t<T>>::value>>
inline decltype(auto) arange(const T & stop) {
  return containers::ArangeContainer<T>(stop);
}

template <class T1, class T2,
          typename = std::enable_if_t<
            std::is_integral<std::common_type_t<T1, T2>>::value>>
inline constexpr decltype(auto) arange(const T1 & start, const T2 & stop) {
  return containers::ArangeContainer<std::common_type_t<T1, T2>>(start, stop);
}

template <class T1, class T2, class T3,
          typename = std::enable_if_t<
            std::is_integral<std::common_type_t<T1, T2, T3>>::value>>
inline constexpr decltype(auto) arange(const T1 & start, const T2 & stop,
                                       const T3 & step) {
  return containers::ArangeContainer<std::common_type_t<T1, T2, T3>>(
      start, stop, step);
}

/* -------------------------------------------------------------------------- */
namespace iterators {
template <class Iterator> class EnumerateIterator {
public:
  using value_type =
      std::tuple<size_t, typename std::iterator_traits<Iterator>::value_type>;
  using difference_type = size_t;
  using pointer =
      std::tuple<size_t, typename std::iterator_traits<Iterator>::pointer>;
  using reference =
      std::tuple<size_t, typename std::iterator_traits<Iterator>::reference>;
  using iterator_category = std::input_iterator_tag;

public:
  explicit EnumerateIterator(Iterator && iterator) : iterator(iterator) {}

  // input iterator ++it
  EnumerateIterator & operator++() {
    ++iterator;
    ++index;
    return *this;
  }

  // input iterator it++
  EnumerateIterator operator++(int) {
    auto cpy = *this;
    this->operator++();
    return cpy;
  }

  // input iterator it != other_it
  bool operator!=(const EnumerateIterator & other) const {
    return iterator != other.iterator;
  }

  // input iterator dereference *it
  decltype(auto) operator*() {
    return std::tuple_cat(std::make_tuple(index), *iterator);
  }

  bool operator==(const EnumerateIterator & other) const {
    return not this->operator!=(other);
  }

private:
  Iterator iterator;
  size_t index{0};
};

template <class Iterator>
inline constexpr decltype(auto) enumerate(Iterator && iterator) {
  return EnumerateIterator<Iterator>(std::forward<Iterator>(iterator));
}

} // namespace iterators

namespace containers {
template <class... Containers> class EnumerateContainer {
public:
  explicit EnumerateContainer(Containers &&... containers)
      : zip_container(std::forward<Containers>(containers)...) {}

  decltype(auto) begin() {
    return iterators::enumerate(zip_container.begin());
  }

  decltype(auto) begin() const {
    return iterators::enumerate(zip_container.begin());
  }

  decltype(auto) end() { return iterators::enumerate(zip_container.end()); }

  decltype(auto) end() const {
    return iterators::enumerate(zip_container.end());
  }

private:
  ZipContainer<Containers...> zip_container;
};
} // namespace containers

template <class... Container>
inline constexpr decltype(auto) enumerate(Container &&... container) {
  return containers::EnumerateContainer<Container...>(
      std::forward<Container>(container)...);
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
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
            std::enable_if_t<details::is_iterator_category_at_least<
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

template <class... zip_container_t>
decltype(auto) make_zip_cat(zip_container_t &&... cont) {
  return make_transform_adaptor(
      zip(std::forward<zip_container_t>(cont)...),
      [](auto && value) { return tuple::flatten(value); });
}

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
        container_begin(container_begin),
        container_it(container_begin) {}

  FilterIterator(const FilterIterator &) = default;

  FilterIterator & operator++() {
    ++filter_it;
    container_it =
        container_begin + *filter_it; // this might return invalid iterators
    return *this;
  }

  decltype(auto) operator*() {
    return *container_it;
  }

  decltype(auto) operator*() const {
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
  container_iterator_t container_it;
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
    static_assert(std::is_same<typename decltype(
                      container.begin())::iterator_category,
                  std::random_access_iterator_tag>::value,
                  "Containers must all have random iterators");
  }

  decltype(auto) begin() const {
    return iterators::make_filter_iterator(filter.begin(),
                                           container.begin());
  }
  decltype(auto) begin() {
    return iterators::make_filter_iterator(filter.begin(),
                                           container.begin());
  }

  decltype(auto) end() const {
    return iterators::make_filter_iterator(filter.end(),
                                           container.begin());
  }
  decltype(auto) end() {
    return iterators::make_filter_iterator(filter.end(),
                                           container.begin());
  }

private:
  filter_t filter;
  Container container;
};
} // namespace containers

template <class filter_t, class Container>
decltype(auto) filter(filter_t && filter,
                      Container && container) {
  return containers::FilterAdaptor<filter_t, Container>(
      std::forward<filter_t>(filter), std::forward<Container>(container));
}

} // namespace AKANTU_ITERATOR_NAMESPACE

namespace std {
template <typename... Its>
struct iterator_traits<::akantu::iterators::ZipIterator<Its...>> {
  using iterator_category = forward_iterator_tag;
  using value_type =
      typename ::akantu::iterators::ZipIterator<Its...>::value_type;
  using difference_type =
      typename ::akantu::iterators::ZipIterator<Its...>::difference_type;
  using pointer = typename ::akantu::iterators::ZipIterator<Its...>::pointer;
  using reference =
      typename ::akantu::iterators::ZipIterator<Its...>::reference;
};

} // namespace std

#endif /* AKANTU_AKA_ITERATORS_HH */