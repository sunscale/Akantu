/**
 * @file   aka_compatibilty_with_cpp_standard.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  Wed Oct 25 2017
 *
 * @brief The content of this file is taken from the possible implementations on
 * http://en.cppreference.com
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
#include <type_traits>
#include <utility>
#include <tuple>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_AKA_COMPATIBILTY_WITH_CPP_STANDARD_HH__
#define __AKANTU_AKA_COMPATIBILTY_WITH_CPP_STANDARD_HH__

namespace aka {

// Part taken from C++14
#if __cplusplus < 201402L
template <bool B, class T = void>
using enable_if_t = typename enable_if<B, T>::type;
#else
template <bool B, class T = void>
using enable_if_t = std::enable_if_t<B, T>;
#endif

// Part taken from C++17
#if __cplusplus < 201703L
// bool_constant
template <bool B> using bool_constant = std::integral_constant<bool, B>;
namespace {
  template <bool B> constexpr bool bool_constant_v = bool_constant<B>::value;
}

// conjunction
template <class...> struct conjunction : std::true_type {};
template <class B1> struct conjunction<B1> : B1 {};
template <class B1, class... Bn>
struct conjunction<B1, Bn...>
    : std::conditional_t<bool(B1::value), conjunction<Bn...>, B1> {};

namespace detail {
  // template <class T>
  // struct is_reference_wrapper : std::false_type {};
  // template <class U>
  // struct is_reference_wrapper<std::reference_wrapper<U>> : std::true_type {};
  // template <class T>
  // constexpr bool is_reference_wrapper_v = is_reference_wrapper<T>::value;

  template <class T, class Type, class T1, class... Args>
  decltype(auto) INVOKE(Type T::*f, T1 && t1, Args &&... args) {
    static_assert(std::is_member_function_pointer<decltype(f)>{} and
                      std::is_base_of<T, std::decay_t<T1>>{},
                  "Does not know what to do with this types");
    return (std::forward<T1>(t1).*f)(std::forward<Args>(args)...);
  }

  // template <class T, class Type, class T1, class... Args>
  // decltype(auto) INVOKE(Type T::*f, T1 && t1, Args &&... args) {
  //   if constexpr (std::is_member_function_pointer_v<decltype(f)>) {
  //     if constexpr (std::is_base_of_v<T, std::decay_t<T1>>)
  //       return (std::forward<T1>(t1).*f)(std::forward<Args>(args)...);
  //     else if constexpr (is_reference_wrapper_v<std::decay_t<T1>>)
  //       return (t1.get().*f)(std::forward<Args>(args)...);
  //     else
  //       return ((*std::forward<T1>(t1)).*f)(std::forward<Args>(args)...);
  //   } else {
  //     static_assert(std::is_member_object_pointer_v<decltype(f)>);
  //     static_assert(sizeof...(args) == 0);
  //     if constexpr (std::is_base_of_v<T, std::decay_t<T1>>)
  //       return std::forward<T1>(t1).*f;
  //     else if constexpr (is_reference_wrapper_v<std::decay_t<T1>>)
  //       return t1.get().*f;
  //     else
  //       return (*std::forward<T1>(t1)).*f;
  //   }
  // }

  template <class F, class... Args>
  decltype(auto) INVOKE(F && f, Args &&... args) {
    return std::forward<F>(f)(std::forward<Args>(args)...);
  }
} // namespace detail

template <class F, class... Args>
decltype(auto) invoke(F && f, Args &&... args) {
  return detail::INVOKE(std::forward<F>(f), std::forward<Args>(args)...);
}

namespace detail {
  template <class F, class Tuple, std::size_t... Is>
  constexpr decltype(auto) apply_impl(F && f, Tuple && t,
                                      std::index_sequence<Is...>) {
    return invoke(std::forward<F>(f),
                  std::get<Is>(std::forward<Tuple>(t))...);
  }
} // namespace detail

template <class F, class Tuple>
constexpr decltype(auto) apply(F && f, Tuple && t) {
  return detail::apply_impl(
      std::forward<F>(f), std::forward<Tuple>(t),
      std::make_index_sequence<std::tuple_size<std::decay_t<Tuple>>::value>{});
}

#endif
} // namespace aka

#endif /* __AKANTU_AKA_COMPATIBILTY_WITH_CPP_STANDARD_HH__ */
