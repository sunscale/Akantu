/**
 * @file   aka_random_generator.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Feb 21 2013
 * @date last modification: Wed Nov 11 2015
 *
 * @brief  generic random generator
 *
 * @section LICENSE
 *
 * Copyright  (©)  2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de Lausanne)
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
#include "aka_array.hh"
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */
#include <ostream>
/* -------------------------------------------------------------------------- */

#if defined(AKANTU_USE_CXX11)
#define __CONST_EXPR constexpr
#else
#define __CONST_EXPR
#endif

#ifndef __AKANTU_AKA_RANDOM_GENERATOR_HH__
#define __AKANTU_AKA_RANDOM_GENERATOR_HH__

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/* List of available distributions                                            */
/* -------------------------------------------------------------------------- */
#define AKANTU_RANDOM_DISTRIBUTION_TYPES                                       \
  ((uniform, UniformDistribution))((weibull, WeibullDistribution))

#define AKANTU_RANDOM_DISTRIBUTION_TYPES_PREFIX(elem) BOOST_PP_CAT(_rdt_, elem)
#define AKANTU_RANDOM_DISTRIBUTION_PREFIX(s, data, elem)                       \
  AKANTU_RANDOM_DISTRIBUTION_TYPES_PREFIX(BOOST_PP_TUPLE_ELEM(2, 0, elem))

enum RandomDistributionType {
  BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(AKANTU_RANDOM_DISTRIBUTION_PREFIX, _,
                                           AKANTU_RANDOM_DISTRIBUTION_TYPES)),
  _rdt_not_defined
};

/* -------------------------------------------------------------------------- */
/* Distribution                                                               */
/* -------------------------------------------------------------------------- */
/// Empty base to be able to store a distribution
template <typename T> class RandomDistributionBase {
public:
  virtual ~RandomDistributionBase() {}
  virtual void printself(std::ostream & stream, int indent = 0) const = 0;
};

/* -------------------------------------------------------------------------- */
/* Uniform distribution                                                       */
/* -------------------------------------------------------------------------- */
template <typename T>
class UniformDistribution : public RandomDistributionBase<T> {
public:
  UniformDistribution(T min = T(0.), T max = T(1.)) : min(min), max(max){};

  /* ------------------------------------------------------------------------ */
  template <template <class> class RandomGenerator>
  T operator()(RandomGenerator<T> & generator) {
    T x = generator() / (RandomGenerator<T>::max() - RandomGenerator<T>::min());
    return (x * (max - min) + min);
  }

  virtual void setParams(std::string value) {
    std::stringstream sstr(value);
    sstr >> min;
    sstr >> max;
  }

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream,
                         __attribute__((unused)) int indent = 0) const {
    stream << "Uniform [ min=" << min << ", max=" << max << " ]";
  };

  /* ------------------------------------------------------------------------ */
private:
  T min;
  T max;
};

/* -------------------------------------------------------------------------- */
/* Weibull distribution                                                       */
/* -------------------------------------------------------------------------- */
template <typename T>
class WeibullDistribution : public RandomDistributionBase<T> {
public:
  WeibullDistribution(T scale, T shape) : m(shape), lambda(scale){};

  /* ------------------------------------------------------------------------ */
  template <template <class> class RandomGenerator>
  T operator()(RandomGenerator<T> & generator) {
    T x = generator() / (RandomGenerator<T>::max() - RandomGenerator<T>::min());
    T e = T(1) / m;
    return lambda * std::pow(-std::log(1. - x), e);
  }

  virtual void setParams(std::string value) {
    std::stringstream sstr(value);
    sstr >> m;
    sstr >> lambda;
  }

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream,
                         __attribute__((unused)) int indent = 0) const {
    stream << "Weibull [ shape=" << m << ", scale=" << lambda << "]";
  }

  /* ------------------------------------------------------------------------ */
private:
  /// shape parameter or Weibull modulus
  T m;
  /// scale parameter
  T lambda;
};

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/* Generator                                                                  */
/* -------------------------------------------------------------------------- */
#if not defined(_WIN32)
template <typename T> class Rand48Generator {
  /* ------------------------------------------------------------------------ */
public:
  inline T operator()() { AKANTU_DEBUG_TO_IMPLEMENT(); }

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const {
    stream << "Rand48Generator [seed=" << _seed << "]";
  }

  /* ------------------------------------------------------------------------ */
public:
  static void seed(long int s) {
    _seed = s;
    srand48(_seed);
  }
  static long int seed() { return _seed; }

  static __CONST_EXPR T min() { return 0.; }
  static __CONST_EXPR T max() { return 1.; }

  /* ------------------------------------------------------------------------ */
private:
  static long int _seed;
};

template <> inline double Rand48Generator<double>::operator()() {
  return drand48();
}
#endif
/* -------------------------------------------------------------------------- */
template <typename T> class RandGenerator {
  /* ------------------------------------------------------------------------ */
public:
  inline T operator()() { return rand(); }

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream,
                         __attribute__((unused)) int indent = 0) const {
    stream << "RandGenerator [seed=" << _seed << "]";
  }

  /* ------------------------------------------------------------------------ */
public:
  static void seed(long int s) {
    _seed = s;
    srand(_seed);
  }
  static long int seed() { return _seed; }

  static __CONST_EXPR T min() { return 0.; }
  static __CONST_EXPR T max() { return RAND_MAX; }

  /* ------------------------------------------------------------------------ */
private:
  static long int _seed;
};

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
template <typename T, template <typename> class Distribution>
class RandomDistributionTypeHelper {
  enum { value = _rdt_not_defined };
};

/* -------------------------------------------------------------------------- */
#define AKANTU_RANDOM_DISTRIBUTION_TYPE_GET_TYPE(r, data, elem)                \
  template <typename T>                                                        \
  struct RandomDistributionTypeHelper<T, BOOST_PP_TUPLE_ELEM(2, 1, elem)> {    \
    enum {                                                                     \
      value = AKANTU_RANDOM_DISTRIBUTION_TYPES_PREFIX(                         \
          BOOST_PP_TUPLE_ELEM(2, 0, elem))                                     \
    };                                                                         \
  };

BOOST_PP_SEQ_FOR_EACH(AKANTU_RANDOM_DISTRIBUTION_TYPE_GET_TYPE, _,
                      AKANTU_RANDOM_DISTRIBUTION_TYPES)

#undef AKANTU_RANDOM_DISTRIBUTION_TYPE_GET_TYPE

/* -------------------------------------------------------------------------- */
/* RandomParameter                                                            */
/* -------------------------------------------------------------------------- */
template <typename T> class RandomParameter {
public:
  RandomParameter(T base_value)
      : base_value(base_value), type(_rdt_not_defined), distribution(NULL) {}

  template <template <typename> class Distribution>
  RandomParameter(T base_value, const Distribution<T> & distribution)
      : base_value(base_value),
        type((RandomDistributionType)
                 RandomDistributionTypeHelper<T, Distribution>::value),
        distribution(new Distribution<T>(distribution)) {}

#define AKANTU_RANDOM_DISTRIBUTION_TYPE_NEW(r, data, elem)                     \
  else if (type == AKANTU_RANDOM_DISTRIBUTION_TYPES_PREFIX(                    \
                       BOOST_PP_TUPLE_ELEM(2, 0, elem))) {                     \
    typedef BOOST_PP_TUPLE_ELEM(2, 1, elem)<T> Dist;                           \
    distribution = new Dist(*static_cast<Dist *>(other.distribution));         \
  }

  RandomParameter(const RandomParameter & other)
      : base_value(other.base_value), type(other.type) {
    if (type == _rdt_not_defined)
      distribution = NULL;
    BOOST_PP_SEQ_FOR_EACH(AKANTU_RANDOM_DISTRIBUTION_TYPE_NEW, _,
                          AKANTU_RANDOM_DISTRIBUTION_TYPES)
  }

  inline void setBaseValue(const T & value) { this->base_value = value; }
  inline T getBaseValue() const { return this->base_value; }

  RandomParameter & operator=(const RandomParameter & other) {
    if (this != &other) {
      base_value = other.base_value;
      type = other.type;
      delete distribution;
      if (type == _rdt_not_defined)
        distribution = NULL;
      BOOST_PP_SEQ_FOR_EACH(AKANTU_RANDOM_DISTRIBUTION_TYPE_NEW, _,
                            AKANTU_RANDOM_DISTRIBUTION_TYPES)
    }
    return *this;
  }
#undef AKANTU_RANDOM_DISTRIBUTION_TYPE_NEW

/* ------------------------------------------------------------------------ */
#define AKANTU_RANDOM_DISTRIBUTION_TYPE_SET(r, data, elem)                     \
  else if (type == AKANTU_RANDOM_DISTRIBUTION_TYPES_PREFIX(                    \
                       BOOST_PP_TUPLE_ELEM(2, 0, elem))) {                     \
    this->set<BOOST_PP_TUPLE_ELEM(2, 1, elem), Generator>(it, end);            \
  }

  template <template <typename> class Generator, class iterator>
  void setValues(iterator it, iterator end) {
    if (type == _rdt_not_defined) {
      for (; it != end; ++it)
        *it = this->base_value;
    }
    BOOST_PP_SEQ_FOR_EACH(AKANTU_RANDOM_DISTRIBUTION_TYPE_SET, _,
                          AKANTU_RANDOM_DISTRIBUTION_TYPES)
  }

#undef AKANTU_RANDOM_DISTRIBUTION_TYPE_SET

  virtual void printself(std::ostream & stream,
                         __attribute__((unused)) int indent = 0) const {
    stream << base_value;
    if (type != _rdt_not_defined)
      stream << " " << *distribution;
  }

private:
  template <template <typename> class Distribution,
            template <typename> class Generator, class iterator>
  void set(iterator it, iterator end) {
    typedef Distribution<T> Dist;
    Dist & dist = *(static_cast<Dist *>(this->distribution));
    Generator<T> gen;
    for (; it != end; ++it)
      *it = this->base_value + dist(gen);
  }

private:
  /// Value with no random variations
  T base_value;

  /// Random distribution type
  RandomDistributionType type;

  /// Random distribution to use
  RandomDistributionBase<T> * distribution;
};

/* -------------------------------------------------------------------------- */
template <typename T>
inline std::ostream & operator<<(std::ostream & stream,
                                 RandomDistributionBase<T> & _this) {
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

__END_AKANTU__

#endif /* __AKANTU_AKA_RANDOM_GENERATOR_HH__ */
