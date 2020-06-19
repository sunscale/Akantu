/**
 * @file   aka_arange_iterator.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  jeu déc 12 2019
 *
 * @brief implementation of arange
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
#include <cstddef>
#include <iterator>
#include <utility>
/* -------------------------------------------------------------------------- */

#ifndef AKA_ARANGE_ITERATOR_HH
#define AKA_ARANGE_ITERATOR_HH

#ifndef AKANTU_ITERATORS_NAMESPACE
#define AKANTU_ITERATORS_NAMESPACE akantu
#endif

namespace AKANTU_ITERATORS_NAMESPACE {

namespace containers {

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

    constexpr T operator[](std::size_t i) {
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
} // namespace AKANTU_ITERATORS_NAMESPACE

namespace std {
template <class T>
struct iterator_traits<
    ::AKANTU_ITERATORS_NAMESPACE::iterators::ArangeIterator<T>> {
private:
  using iterator_type =
      typename ::AKANTU_ITERATORS_NAMESPACE::iterators::ArangeIterator<T>;

public:
  using iterator_category = typename iterator_type::iterator_category;
  using value_type = typename iterator_type::value_type;
  using difference_type = typename iterator_type::difference_type;
  using pointer = typename iterator_type::pointer;
  using reference = typename iterator_type::reference;
};

} // namespace std

#endif /* AKA_ARANGE_ITERATOR_HH */
