/**
 * @file   aka_random_generator.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Feb 21 2013
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  generic random generator
 *
 * @section LICENSE
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_array.hh"
/* -------------------------------------------------------------------------- */
#include <random>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_AKA_RANDOM_GENERATOR_HH__
#define __AKANTU_AKA_RANDOM_GENERATOR_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
/* List of available distributions                                            */
/* -------------------------------------------------------------------------- */
// clang-format off
#define AKANTU_RANDOM_DISTRIBUTION_TYPES                \
  ((uniform      , std::uniform_real_distribution ))    \
  ((exponential  , std::exponential_distribution  ))    \
  ((gamma        , std::gamma_distribution        ))    \
  ((weibull      , std::weibull_distribution      ))    \
  ((extreme_value, std::extreme_value_distribution))    \
  ((normal       , std::normal_distribution       ))    \
  ((lognormal    , std::lognormal_distribution    ))    \
  ((chi_squared  , std::chi_squared_distribution  ))    \
  ((cauchy       , std::cauchy_distribution       ))    \
  ((fisher_f     , std::fisher_f_distribution     ))    \
  ((student_t    , std::student_t_distribution    ))
// clang-format on

#define AKANTU_RANDOM_DISTRIBUTION_TYPES_PREFIX(elem) BOOST_PP_CAT(_rdt_, elem)
#define AKANTU_RANDOM_DISTRIBUTION_PREFIX(s, data, elem)                       \
  AKANTU_RANDOM_DISTRIBUTION_TYPES_PREFIX(BOOST_PP_TUPLE_ELEM(2, 0, elem))

enum RandomDistributionType {
  BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(AKANTU_RANDOM_DISTRIBUTION_PREFIX, _,
                                           AKANTU_RANDOM_DISTRIBUTION_TYPES)),
  _rdt_not_defined
};

/* -------------------------------------------------------------------------- */
/* Generator                                                                  */
/* -------------------------------------------------------------------------- */
template <typename T> class RandomGenerator {
  /* ------------------------------------------------------------------------ */
private:
  static long int _seed;
  static std::default_random_engine generator;
  /* ------------------------------------------------------------------------ */
public:
  inline T operator()() { return generator(); }

  /// function to print the contain of the class
  void printself(std::ostream & stream, int) const {
    stream << "RandGenerator [seed=" << _seed << "]";
  }

  /* ------------------------------------------------------------------------ */
public:
  static void seed(long int s) {
    _seed = s;
    generator.seed(_seed);
  }
  static long int seed() { return _seed; }

  static constexpr T min() { return generator.min(); }
  static constexpr T max() { return generator.max(); }
};

template <typename T> long int RandomGenerator<T>::_seed;
template <typename T> std::default_random_engine RandomGenerator<T>::generator;

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
#undef AKANTU_RANDOM_DISTRIBUTION_PREFIX

#define AKANTU_RANDOM_DISTRIBUTION_TYPE_PRINT_CASE(r, data, elem)              \
  case AKANTU_RANDOM_DISTRIBUTION_TYPES_PREFIX(                                \
      BOOST_PP_TUPLE_ELEM(2, 0, elem)): {                                      \
    stream << BOOST_PP_STRINGIZE(AKANTU_RANDOM_DISTRIBUTION_TYPES_PREFIX(      \
        BOOST_PP_TUPLE_ELEM(2, 0, elem)));                                     \
    break;                                                                     \
  }

inline std::ostream & operator<<(std::ostream & stream,
                                 RandomDistributionType type) {
  switch (type) {
    BOOST_PP_SEQ_FOR_EACH(AKANTU_RANDOM_DISTRIBUTION_TYPE_PRINT_CASE, _,
                          AKANTU_RANDOM_DISTRIBUTION_TYPES)
  default:
    stream << UInt(type) << " not a RandomDistributionType";
    break;
  }
  return stream;
}
#undef AKANTU_RANDOM_DISTRIBUTION_TYPE_PRINT_CASE

/* -------------------------------------------------------------------------- */
/* Some Helper                                                                */
/* -------------------------------------------------------------------------- */
template <typename T, class Distribution> class RandomDistributionTypeHelper {
  enum { value = _rdt_not_defined };
};

/* -------------------------------------------------------------------------- */
#define AKANTU_RANDOM_DISTRIBUTION_TYPE_GET_TYPE(r, data, elem)                \
  template <typename T>                                                        \
      struct RandomDistributionTypeHelper<T, BOOST_PP_TUPLE_ELEM(2, 1, elem) < \
                                                 T>> {                         \
    enum {                                                                     \
      value = AKANTU_RANDOM_DISTRIBUTION_TYPES_PREFIX(                         \
          BOOST_PP_TUPLE_ELEM(2, 0, elem))                                     \
    };                                                                         \
                                                                               \
    static void printself(std::ostream & stream) {                             \
      stream << BOOST_PP_STRINGIZE(BOOST_PP_TUPLE_ELEM(2, 0, elem));           \
    }                                                                          \
  };

BOOST_PP_SEQ_FOR_EACH(AKANTU_RANDOM_DISTRIBUTION_TYPE_GET_TYPE, _,
                      AKANTU_RANDOM_DISTRIBUTION_TYPES)

#undef AKANTU_RANDOM_DISTRIBUTION_TYPE_GET_TYPE

/* -------------------------------------------------------------------------- */
template <class T> class RandomDistribution {
public:
  virtual ~RandomDistribution() = default;
  virtual T operator()(RandomGenerator<UInt> & gen) = 0;
  virtual std::unique_ptr<RandomDistribution<T>> make_unique() const = 0;
  virtual void printself(std::ostream & stream, int = 0) const = 0;
};

template <class T, class Distribution>
class RandomDistributionProxy : public RandomDistribution<T> {
public:
  explicit RandomDistributionProxy(Distribution dist)
      : distribution(std::move(dist)) {}

  T operator()(RandomGenerator<UInt> & gen) override {
    return distribution(gen);
  }

  std::unique_ptr<RandomDistribution<T>> make_unique() const override {
    return std::make_unique<RandomDistributionProxy<T, Distribution>>(
        distribution);
  }

  void printself(std::ostream & stream, int = 0) const override {
    RandomDistributionTypeHelper<T, Distribution>::printself(stream);
    stream << " [ " << distribution << " ]";
  }

private:
  Distribution distribution;
};

/* -------------------------------------------------------------------------- */
/* RandomParameter                                                            */
/* -------------------------------------------------------------------------- */
template <typename T> class RandomParameter {
public:
  template <class Distribution>
  explicit RandomParameter(T base_value, Distribution dist)
      : base_value(base_value),
        type(RandomDistributionType(
            RandomDistributionTypeHelper<T, Distribution>::value)),
        distribution_proxy(
            std::make_unique<RandomDistributionProxy<T, Distribution>>(
                std::move(dist))) {}

  explicit RandomParameter(T base_value)
      : base_value(base_value),
        type(RandomDistributionType(
            RandomDistributionTypeHelper<
                T, std::uniform_real_distribution<T>>::value)),
        distribution_proxy(
            std::make_unique<
                RandomDistributionProxy<T, std::uniform_real_distribution<T>>>(
                std::uniform_real_distribution<T>(0., 0.))) {}

  RandomParameter(const RandomParameter & other)
      : base_value(other.base_value), type(other.type),
        distribution_proxy(other.distribution_proxy->make_unique()) {}

  RandomParameter & operator=(const RandomParameter & other) {
    distribution_proxy = other.distribution_proxy->make_unique();
    base_value = other.base_value;
    type = other.type;
    return *this;
  }

  virtual ~RandomParameter() = default;

  inline void setBaseValue(const T & value) { this->base_value = value; }
  inline T getBaseValue() const { return this->base_value; }

  template <template <typename> class Generator, class iterator>
  void setValues(iterator it, iterator end) {
    RandomGenerator<UInt> gen;
    for (; it != end; ++it)
      *it = this->base_value + (*distribution_proxy)(gen);
  }

  virtual void printself(std::ostream & stream,
                         __attribute__((unused)) int indent = 0) const {
    stream << base_value;
    stream << " + " << *distribution_proxy;
  }

private:
  /// Value with no random variations
  T base_value;

  /// Random distribution type
  RandomDistributionType type;

  /// Proxy to store a std random distribution
  std::unique_ptr<RandomDistribution<T>> distribution_proxy;
};

/* -------------------------------------------------------------------------- */
template <typename T>
inline std::ostream & operator<<(std::ostream & stream,
                                 RandomDistribution<T> & _this) {
  _this.printself(stream);
  return stream;
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline std::ostream & operator<<(std::ostream & stream,
                                 RandomParameter<T> & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#endif /* __AKANTU_AKA_RANDOM_GENERATOR_HH__ */
