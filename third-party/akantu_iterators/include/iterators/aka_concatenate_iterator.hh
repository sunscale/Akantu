/**
 * @file   aka_concatenate_iterator.hh
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
#include "aka_compatibilty_with_cpp_standard.hh"
#include "aka_iterator_tools.hh"
#include "aka_tuple_tools.hh"
/* -------------------------------------------------------------------------- */
#include <iterator>
#include <utility>
/* -------------------------------------------------------------------------- */

#ifndef AKA_CONCATENATE_ITERATOR_H
#define AKA_CONCATENATE_ITERATOR_H

#ifndef AKANTU_ITERATORS_NAMESPACE
#define AKANTU_ITERATORS_NAMESPACE akantu
#endif

namespace AKANTU_ITERATORS_NAMESPACE {

namespace iterators {

  /* ------------------------------------------------------------------------ */
  template <class... Iterators>
  class ConcatIterator
      : public details::CopyAssignmentEnabler<
            aka::conjunction<std::is_copy_assignable<Iterators>...,
                             std::is_copy_constructible<Iterators>...>::value>,
        public details::MoveAssignmentEnabler<
            aka::conjunction<std::is_move_assignable<Iterators>...,
                             std::is_move_constructible<Iterators>...>::value> {
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
    using iterator_category = std::input_iterator_tag;

  public:
    explicit ConcatIterator(tuple_t iterators, tuple_t end_iterators)
        : iterators(std::move(iterators)),
          end_iterators(std::move(end_iterators)) {}

    /* ---------------------------------------------------------------------- */
    // input iterator ++it
    ConcatIterator & operator++() {
      auto && ends =
          tuple::transform([](auto && a, auto && b) { return a == b; },
                           iterators, end_iterators);
      auto && pos = tuple::find(ends, false);
      ++(tuple::dynamic_get(pos, iterators));
      return *this;
    }

    // input iterator it++
    ConcatIterator operator++(int) {
      auto cpy = *this;
      this->operator++();
      return cpy;
    }

    // input iterator it != other_it
    bool operator!=(const ConcatIterator & other) const {
      return tuple::are_not_equal(iterators, other.iterators);
    }

    // input iterator dereference *it
    decltype(auto) operator*() {
      auto && ends =
          tuple::transform([](auto && a, auto && b) { return a == b; },
                           iterators, end_iterators);
      auto && pos = tuple::find(ends, false);
      return *(tuple::dynamic_get(pos, iterators));
    }

    template <
        class iterator_category_ = iterator_category,
        std::enable_if_t<aka::is_iterator_category_at_least<
            iterator_category_, std::forward_iterator_tag>::value> * = nullptr>
    bool operator==(const ConcatIterator & other) const {
      return not tuple::are_not_equal(iterators, other.iterators);
    }

  private:
    tuple_t iterators;
    tuple_t end_iterators;
  };

} // namespace iterators

/* -------------------------------------------------------------------------- */
template <class... Iterators>
decltype(auto)
concat_iterator(std::tuple<Iterators...> && iterators_tuple,
                std::tuple<Iterators...> && end_iterators_tuple) {
  auto concat = iterators::ConcatIterator<Iterators...>(
      std::forward<decltype(iterators_tuple)>(iterators_tuple),
      std::forward<decltype(end_iterators_tuple)>(end_iterators_tuple));
  return concat;
}

/* -------------------------------------------------------------------------- */
namespace containers {
  template <class... Containers> class ConcatContainer {
    using containers_t = std::tuple<Containers...>;

  public:
    explicit ConcatContainer(Containers &&... containers)
        : containers(std::forward<Containers>(containers)...) {}

    decltype(auto) begin() const {
      return concat_iterator(
          tuple::transform([](auto && c) { return c.begin(); },
                           std::forward<containers_t>(containers)),
          tuple::transform([](auto && c) { return c.end(); },
                           std::forward<containers_t>(containers)));
    }

    decltype(auto) end() const {
      return concat_iterator(
          tuple::transform([](auto && c) { return c.end(); },
                           std::forward<containers_t>(containers)),
          tuple::transform([](auto && c) { return c.end(); },
                           std::forward<containers_t>(containers)));
    }

    decltype(auto) begin() {
       return concat_iterator(
          tuple::transform([](auto && c) { return c.begin(); },
                           std::forward<containers_t>(containers)),
          tuple::transform([](auto && c) { return c.end(); },
                           std::forward<containers_t>(containers)));
    }

    decltype(auto) end() {
      return concat_iterator(
          tuple::transform([](auto && c) { return c.end(); },
                           std::forward<containers_t>(containers)),
          tuple::transform([](auto && c) { return c.end(); },
                           std::forward<containers_t>(containers)));
    }

  private:
    containers_t containers;
  };

  /* ------------------------------------------------------------------------ */
  template <class... Containers> decltype(auto) concat(Containers &&... conts) {
    return containers::ConcatContainer<Containers...>(
        std::forward<Containers>(conts)...);
  }

} // namespace containers
} // namespace AKANTU_ITERATORS_NAMESPACE

namespace std {
template <typename... Its>
struct iterator_traits<
    ::AKANTU_ITERATORS_NAMESPACE::iterators::ConcatIterator<Its...>> {
private:
  using iterator_type =
      typename ::AKANTU_ITERATORS_NAMESPACE::iterators::ConcatIterator<Its...>;

public:
  using iterator_category = typename iterator_type::iterator_category;
  using value_type = typename iterator_type::value_type;
  using difference_type = typename iterator_type::difference_type;
  using pointer = typename iterator_type::pointer;
  using reference = typename iterator_type::reference;
};
} // namespace std

#endif // __AKA_CONCATENATE_ITERATOR_H_
