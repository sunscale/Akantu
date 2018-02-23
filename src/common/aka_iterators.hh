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
#include <tuple>
#include <utility>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_AKA_ITERATORS_HH__
#define __AKANTU_AKA_ITERATORS_HH__

namespace akantu {

namespace tuple {
  /* ------------------------------------------------------------------------ */
  namespace details {
    template <size_t N> struct Foreach {
      template <class Tuple>
      static inline bool not_equal(Tuple && a, Tuple && b) {
        if (std::get<N - 1>(std::forward<Tuple>(a)) ==
            std::get<N - 1>(std::forward<Tuple>(b)))
          return false;
        return Foreach<N - 1>::not_equal(std::forward<Tuple>(a),
                                         std::forward<Tuple>(b));
      }
    };

    /* ---------------------------------------------------------------------- */
    template <> struct Foreach<0> {
      template <class Tuple>
      static inline bool not_equal(Tuple && a, Tuple && b) {
        return std::get<0>(std::forward<Tuple>(a)) !=
               std::get<0>(std::forward<Tuple>(b));
      }
    };

    template <typename... Ts>
    decltype(auto) make_tuple_no_decay(Ts &&... args) {
      return std::tuple<Ts...>(std::forward<Ts>(args)...);
    }

    template <class F, class Tuple, size_t... Is>
    void foreach_impl(F && func, Tuple && tuple,
                      std::index_sequence<Is...> &&) {
      (void)std::initializer_list<int>{
          (std::forward<F>(func)(std::get<Is>(std::forward<Tuple>(tuple))),
           0)...};
    }

    template <class F, class Tuple, size_t... Is>
    decltype(auto) transform_impl(F && func, Tuple && tuple,
                                  std::index_sequence<Is...> &&) {
      return make_tuple_no_decay(
          std::forward<F>(func)(std::get<Is>(std::forward<Tuple>(tuple)))...);
    }
  } // namespace details

  /* ------------------------------------------------------------------------ */
  template <class Tuple> bool are_not_equal(Tuple && a, Tuple && b) {
    return details::Foreach<std::tuple_size<std::decay_t<Tuple>>::value>::
        not_equal(std::forward<Tuple>(a), std::forward<Tuple>(b));
  }

  template <class F, class Tuple> void foreach (F && func, Tuple && tuple) {
    return details::foreach_impl(
        std::forward<F>(func), std::forward<Tuple>(tuple),
        std::make_index_sequence<
            std::tuple_size<std::decay_t<Tuple>>::value>{});
  }

  template <class F, class Tuple>
  decltype(auto) transform(F && func, Tuple && tuple) {
    return details::transform_impl(
        std::forward<F>(func), std::forward<Tuple>(tuple),
        std::make_index_sequence<
            std::tuple_size<std::decay_t<Tuple>>::value>{});
  }
} // namespace tuple

/* -------------------------------------------------------------------------- */
namespace iterators {
  namespace {
    template <typename It, typename category>
    struct is_at_least_category
        : std::is_same<std::common_type_t<
                           typename std::iterator_traits<It>::iterator_category,
                           category>,
                       category> {};
  }

  template <class... Iterators> class ZipIterator {
  public:
    using value_type =
        std::tuple<typename std::iterator_traits<Iterators>::value_type...>;
    using difference_type = std::common_type_t<
        typename std::iterator_traits<Iterators>::difference_type...>;
    using pointer =
        std::tuple<typename std::iterator_traits<Iterators>::pointer...>;
    using reference =
        std::tuple<typename std::iterator_traits<Iterators>::reference...>;
    using iterator_category = std::input_iterator_tag;
    // std::common_type_t<typename
    // std::iterator_traits<Iterators>::iterator_category...>;

  private:
    using tuple_t = std::tuple<Iterators...>;

  public:
    explicit ZipIterator(tuple_t iterators) : iterators(std::move(iterators)) {}

    // input iterator ++it
    ZipIterator & operator++() {
      tuple::foreach ([](auto && it) { ++it; }, iterators);
      return *this;
    }

    // input iterator it++
    ZipIterator operator++(int) {
      auto cpy = *this;
      tuple::foreach ([](auto && it) { ++it; }, iterators);
      return cpy;
    }

    // input iterator it != other_it
    bool operator!=(const ZipIterator & other) const {
      return tuple::are_not_equal(iterators, other.iterators);
    }

    // input iterator dereference *it
    decltype(auto) operator*() {
      return tuple::transform([](auto && it) -> decltype(auto) { return *it; },
                              iterators);
    }

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

  private:
    containers_t containers;
  };

  template <class Iterator> class Range {
  public:
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
    using iterator_category = std::input_iterator_tag;

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
template <class Container>
inline constexpr decltype(auto) enumerate(Container && container,
                                          size_t start_ = 0) {
  auto stop = std::forward<Container>(container).size();
  decltype(stop) start = start_;
  return zip(arange(start, stop), std::forward<Container>(container));
}

} // namespace akantu

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
}

#endif /* __AKANTU_AKA_ITERATORS_HH__ */
