/**
 * @file   test_gtest_utils.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Nov 14 2017
 * @date last modification: Wed Feb 21 2018
 *
 * @brief  Utils to help write tests
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
#include "aka_common.hh"
#include "aka_iterators.hh"
/* -------------------------------------------------------------------------- */
#include <boost/preprocessor.hpp>
#include <gtest/gtest.h>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_TEST_GTEST_UTILS_HH__
#define __AKANTU_TEST_GTEST_UTILS_HH__

#if !defined(TYPED_TEST_SUITE)
#define TYPED_TEST_SUITE(...) TYPED_TEST_CASE(__VA_ARGS__)
#endif

#if !defined(TYPED_TEST_SUITE_P)
#define TYPED_TEST_SUITE_P(...) TYPED_TEST_CASE_P(__VA_ARGS__)
#endif

#if !defined(REGISTER_TYPED_TEST_SUITE_P)
#define REGISTER_TYPED_TEST_SUITE_P(...) REGISTER_TYPED_TEST_CASE_P(__VA_ARGS__)
#endif

#if !defined(INSTANTIATE_TYPED_TEST_SUITE_P)
#define INSTANTIATE_TYPED_TEST_SUITE_P(...)                                    \
  INSTANTIATE_TYPED_TEST_CASE_P(__VA_ARGS__)
#endif

namespace {

/* -------------------------------------------------------------------------- */
template <::akantu::ElementType t>
using element_type_t = std::integral_constant<::akantu::ElementType, t>;

/* -------------------------------------------------------------------------- */
template <typename... T> struct gtest_list {};

template <typename... Ts> struct gtest_list<std::tuple<Ts...>> {
  using type = ::testing::Types<Ts...>;
};

template <typename... T> using gtest_list_t = typename gtest_list<T...>::type;

/* -------------------------------------------------------------------------- */
template <typename... T> struct tuple_concat {};

template <typename... T1s, typename... T2s>
struct tuple_concat<std::tuple<T1s...>, std::tuple<T2s...>> {
  using type = std::tuple<T1s..., T2s...>;
};

template <typename... T>
using tuple_concat_t = typename tuple_concat<T...>::type;

/* -------------------------------------------------------------------------- */
template <template <typename> class Pred, typename... Ts>
struct tuple_filter {};

template <template <typename> class Pred, typename T>
struct tuple_filter<Pred, std::tuple<T>> {
  using type = std::conditional_t<Pred<T>::value, std::tuple<T>, std::tuple<>>;
};

template <template <typename> class Pred, typename T, typename... Ts>
struct tuple_filter<Pred, std::tuple<T, Ts...>> {
  using type =
      tuple_concat_t<typename tuple_filter<Pred, std::tuple<T>>::type,
                     typename tuple_filter<Pred, std::tuple<Ts...>>::type>;
};

template <template <typename> class Pred, typename... Ts>
using tuple_filter_t = typename tuple_filter<Pred, Ts...>::type;

/* -------------------------------------------------------------------------- */
template <size_t N, typename... Ts> struct tuple_split {};

template <size_t N, typename T, typename... Ts>
struct tuple_split<N, std::tuple<T, Ts...>> {
protected:
  using split = tuple_split<N - 1, std::tuple<Ts...>>;

public:
  using type = tuple_concat_t<std::tuple<T>, typename split::type>;
  using type_tail = typename split::type_tail;
};

template <typename T, typename... Ts>
struct tuple_split<1, std::tuple<T, Ts...>> {
  using type = std::tuple<T>;
  using type_tail = std::tuple<Ts...>;
};

template <size_t N, typename... T>
using tuple_split_t = typename tuple_split<N, T...>::type;

template <size_t N, typename... T>
using tuple_split_tail_t = typename tuple_split<N, T...>::type_tail;

/* -------------------------------------------------------------------------- */
template <typename... T> struct cross_product {};

template <typename... T2s>
struct cross_product<std::tuple<>, std::tuple<T2s...>> {
  using type = std::tuple<>;
};

template <typename T1, typename... T1s, typename... T2s>
struct cross_product<std::tuple<T1, T1s...>, std::tuple<T2s...>> {
  using type = tuple_concat_t<
      std::tuple<std::tuple<T1, T2s>...>,
      typename cross_product<std::tuple<T1s...>, std::tuple<T2s...>>::type>;
};

template <typename... T>
using cross_product_t = typename cross_product<T...>::type;
/* -------------------------------------------------------------------------- */

} // namespace

#define OP_CAT(s, data, elem) BOOST_PP_CAT(_element_type, elem)

// creating a type instead of a using helps to debug
#define AKANTU_DECLARE_ELEMENT_TYPE_STRUCT(r, data, elem)                      \
  struct BOOST_PP_CAT(_element_type, elem)                                     \
      : public element_type_t<::akantu::elem> {};

BOOST_PP_SEQ_FOR_EACH(AKANTU_DECLARE_ELEMENT_TYPE_STRUCT, _,
                      AKANTU_ALL_ELEMENT_TYPE)

#undef AKANTU_DECLARE_ELEMENT_TYPE_STRUCT

using TestElementTypesAll = std::tuple<BOOST_PP_SEQ_ENUM(
    BOOST_PP_SEQ_TRANSFORM(OP_CAT, _, AKANTU_ek_regular_ELEMENT_TYPE))>;

#if defined(AKANTU_COHESIVE_ELEMENT)
using TestCohesiveElementTypes = std::tuple<BOOST_PP_SEQ_ENUM(
    BOOST_PP_SEQ_TRANSFORM(OP_CAT, _, AKANTU_ek_cohesive_ELEMENT_TYPE))>;
#endif

#if defined(AKANTU_STRUCTURAL_MECHANICS)
using TestElementTypesStructural = std::tuple<BOOST_PP_SEQ_ENUM(
    BOOST_PP_SEQ_TRANSFORM(OP_CAT, _, AKANTU_ek_structural_ELEMENT_TYPE))>;
#endif

using TestAllDimensions = std::tuple<std::integral_constant<unsigned int, 1>,
                                     std::integral_constant<unsigned int, 2>,
                                     std::integral_constant<unsigned int, 3>>;

template <typename T, ::akantu::ElementType type>
using is_element = aka::bool_constant<T::value == type>;

template <typename T>
using not_is_point_1 = aka::negation<is_element<T, ::akantu::_point_1>>;

using TestElementTypes = tuple_filter_t<not_is_point_1, TestElementTypesAll>;

#if defined(AKANTU_STRUCTURAL_MECHANICS)
using StructuralTestElementTypes =
    tuple_filter_t<not_is_point_1, TestElementTypesStructural>;
#endif

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
template <size_t degree> class Polynomial {
public:
  Polynomial() = default;

  Polynomial(std::initializer_list<double> && init) {
    for (auto && pair : akantu::zip(init, constants))
      std::get<1>(pair) = std::get<0>(pair);
  }

  double operator()(double x) {
    double res = 0.;
    for (auto && vals : akantu::enumerate(constants)) {
      double a;
      int k;
      std::tie(k, a) = vals;
      res += a * std::pow(x, k);
    }
    return res;
  }

  Polynomial extract(size_t pdegree) {
    Polynomial<degree> extract(*this);
    for (size_t d = pdegree + 1; d < degree + 1; ++d)
      extract.constants[d] = 0;
    return extract;
  }

  auto integral() {
    Polynomial<degree + 1> integral_;
    integral_.set(0, 0.);
    ;
    for (size_t d = 0; d < degree + 1; ++d) {
      integral_.set(1 + d, get(d) / double(d + 1));
    }
    return integral_;
  }

  auto integrate(double a, double b) {
    auto primitive = integral();
    return (primitive(b) - primitive(a));
  }

  double get(int i) const { return constants[i]; }

  void set(int i, double a) { constants[i] = a; }

protected:
  std::array<double, degree + 1> constants;
};

template <size_t degree>
std::ostream & operator<<(std::ostream & stream, const Polynomial<degree> & p) {
  for (size_t d = 0; d < degree + 1; ++d) {
    if (d != 0)
      stream << " + ";

    stream << p.get(degree - d);
    if (d != degree)
      stream << "x ^ " << degree - d;
  }
  return stream;
}

/* -------------------------------------------------------------------------- */

#endif /* __AKANTU_TEST_GTEST_UTILS_HH__ */
