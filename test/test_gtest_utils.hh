/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_iterators.hh"
/* -------------------------------------------------------------------------- */
#include <boost/preprocessor.hpp>
#include <gtest/gtest.h>
#include <tuple>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_TEST_GTEST_UTILS_HH__
#define __AKANTU_TEST_GTEST_UTILS_HH__

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
template <template <typename> class Pred, typename... Ts> struct tuple_filter {};

template <template <typename> class Pred, typename T>
struct tuple_filter<Pred, std::tuple<T>> {
  using type = std::conditional_t<Pred<T>::value, std::tuple<T>, std::tuple<>>;
};

template <template <typename> class Pred, typename T, typename... Ts>
struct tuple_filter<Pred, std::tuple<T, Ts...>> {
  using type = tuple_concat_t<
    typename tuple_filter<Pred, std::tuple<T>>::type,
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

#define OP_CAT(s, data, elem) element_type_t<::akantu::elem>

using TestElementTypesAll = std::tuple<BOOST_PP_SEQ_ENUM(
    BOOST_PP_SEQ_TRANSFORM(OP_CAT, _, AKANTU_ek_regular_ELEMENT_TYPE))>;
} // namespace

template <typename T>
using is_point_1 = std::is_same<T, element_type_t<::akantu::_point_1>>;

template <typename T>
using not_is_point_1 =
    aka::negation<std::is_same<T, element_type_t<::akantu::_point_1>>>;

using TestElementTypes = tuple_filter_t<not_is_point_1, TestElementTypesAll>;

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
