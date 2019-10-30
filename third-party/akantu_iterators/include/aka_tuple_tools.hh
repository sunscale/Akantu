/**
 * @file   aka_tuple_tools.hh
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
 * Copyright 2019 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_AKA_TUPLE_TOOLS_HH
#define AKANTU_AKA_TUPLE_TOOLS_HH

#ifndef AKANTU_ITERATOR_NAMESPACE
#define AKANTU_ITERATOR_NAMESPACE akantu
#endif

namespace AKANTU_ITERATOR_NAMESPACE {

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

  namespace details {
    template <class Tuple, std::size_t... Is>
    decltype(auto) flatten(Tuple && tuples, std::index_sequence<Is...>) {
      return std::tuple_cat(std::get<Is>(tuples)...);
    }
  } // namespace details

  template <class Tuple> decltype(auto) flatten(Tuple && tuples) {
    return details::flatten(std::forward<Tuple>(tuples),
                            std::make_index_sequence<
                                std::tuple_size<std::decay_t<Tuple>>::value>());
  }

} // namespace tuple

} // namespace AKANTU_ITERATOR_NAMESPACE

#endif /* AKANTU_AKA_TUPLE_TOOLS_HH */
