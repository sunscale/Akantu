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
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_AKA_COMPATIBILTY_WITH_CPP_STANDARD_HH__
#define __AKANTU_AKA_COMPATIBILTY_WITH_CPP_STANDARD_HH__

namespace std {

// Part taken from C++14
#if __cplusplus < 201402L
template <bool B, class T = void>
using enable_if_t = typename enable_if<B, T>::type;
#endif

// Part taken from C++17
#if __cplusplus < 201703L
// bool_constant
template <bool B> using bool_constant = integral_constant<bool, B>;
namespace {
  template <bool B> constexpr bool bool_constant_v = bool_constant<B>::value;
}

// conjunction
template <class...> struct conjunction : std::true_type {};
template <class B1> struct conjunction<B1> : B1 {};
template <class B1, class... Bn>
struct conjunction<B1, Bn...>
    : std::conditional_t<bool(B1::value), conjunction<Bn...>, B1> {};

#endif
}

#endif /* __AKANTU_AKA_COMPATIBILTY_WITH_CPP_STANDARD_HH__ */
