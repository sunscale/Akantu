/**
 * @file   algebraic_parser.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Nov 13 2013
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  algebraic_parser definition of the grammar
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
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */
// Boost
#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/include/qi.hpp>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_ALGEBRAIC_PARSER_HH__
#define __AKANTU_ALGEBRAIC_PARSER_HH__

namespace spirit = boost::spirit;
namespace qi = boost::spirit::qi;
namespace lbs = boost::spirit::qi::labels;
namespace ascii = boost::spirit::ascii;
namespace phx = boost::phoenix;

namespace akantu {
namespace parser {
  struct algebraic_error_handler_ {
    template <typename, typename, typename> struct result {
      using type = void;
    };

    template <typename Iterator>
    void operator()(qi::info const & what, Iterator err_pos,
                    Iterator last) const {
      AKANTU_EXCEPTION(
          "Error! Expecting "
          << what // what failed?
          << " here: \""
          << std::string(err_pos, last) // iterators to error-pos, end
          << "\"");
    }
  };

  static Real my_min(Real a, Real b) { return std::min(a, b); }
  static Real my_max(Real a, Real b) { return std::max(a, b); }
  static Real my_pow(Real a, Real b) { return std::pow(a, b); }

  static Real eval_param(const ID & a, const ParserSection & section) {
    return section.getParameter(a, _ppsc_current_and_parent_scope);
  }

  static Real unary_func(Real (*func)(Real), Real a) { return func(a); }

  static Real binary_func(Real (*func)(Real, Real), Real a, Real b) {
    return func(a, b);
  }

  template <class Iterator, typename Skipper = spirit::unused_type>
  struct AlgebraicGrammar : qi::grammar<Iterator, Real(), Skipper> {
    AlgebraicGrammar(const ParserSection & section)
        : AlgebraicGrammar::base_type(start, "algebraic_grammar"),
          section(section) {
      // phx::function<lazy_pow_>  lazy_pow;
      // phx::function<lazy_unary_func_>  lazy_unary_func;
      // phx::function<lazy_binary_func_> lazy_binary_func;
      // phx::function<lazy_eval_param_> lazy_eval_param;
      /* clang-format off */
      start
        =   expr.alias()
        ;

      expr
        =   term                  [ lbs::_val  = lbs::_1 ]
            >> *( ('+' > term     [ lbs::_val += lbs::_1 ])
                | ('-' > term     [ lbs::_val -= lbs::_1 ])
                )
        ;

      term
        =   factor                [ lbs::_val  = lbs::_1 ]
            >> *( ('*' > factor   [ lbs::_val *= lbs::_1 ])
                | ('/' > factor   [ lbs::_val /= lbs::_1 ])
                )
        ;

      factor
        =   number                [ lbs::_val = lbs::_1 ]
        >> *("**" > number    [ lbs::_val = phx::bind(&my_pow, lbs::_val, lbs::_1) ])
        ;

      number
        =   real                  [ lbs::_val =  lbs::_1 ]
        |   ('-' > number         [ lbs::_val = -lbs::_1 ])
        |   ('+' > number         [ lbs::_val =  lbs::_1 ])
        |   constant              [ lbs::_val =  lbs::_1 ]
        |   function              [ lbs::_val =  lbs::_1 ]
        |   ('(' > expr > ')')    [ lbs::_val =  lbs::_1 ]
        |   variable              [ lbs::_val =  lbs::_1 ]
        ;

      function
        =   (qi::no_case[unary_function]
             > '('
             > expr
             > ')')               [ lbs::_val = phx::bind(&unary_func, lbs::_1, lbs::_2) ]
        |   (qi::no_case[binary_function]
             > '(' >> expr
             > ',' >> expr
             > ')')               [ lbs::_val = phx::bind(&binary_func ,lbs::_1, lbs::_2, lbs::_3) ]
        ;

      variable
        =   key [ lbs::_val = phx::bind(&eval_param, lbs::_1, section) ]
        ;

      key
        =   qi::no_skip[qi::char_("a-zA-Z_") >> *qi::char_("a-zA-Z_0-9")] // coming from the InputFileGrammar
        ;

#ifndef M_PI
#  define M_PI          3.14159265358979323846
#endif
#ifndef M_E
#  define M_E           2.7182818284590452354
#endif
      constant.add
        ("pi", M_PI)
        ("e",  M_E);

      unary_function.add
        ("abs"   , &std::abs   )
        ("acos"  , &std::acos  )
        ("asin"  , &std::asin  )
        ("atan"  , &std::atan  )
        ("ceil"  , &std::ceil  )
        ("cos"   , &std::cos   )
        ("cosh"  , &std::cosh  )
        ("exp"   , &std::exp   )
        ("floor" , &std::floor )
        ("log10" , &std::log10 )
        ("log"   , &std::log   )
        ("sin"   , &std::sin   )
        ("sinh"  , &std::sinh  )
        ("sqrt"  , &std::sqrt  )
        ("tan"   , &std::tan   )
        ("tanh"  , &std::tanh  )
        ("acosh" , &std::acosh )
        ("asinh" , &std::asinh )
        ("atanh" , &std::atanh )
        ("exp2"  , &std::exp2  )
        ("expm1" , &std::expm1 )
        ("log1p" , &std::log1p )
        ("log2"  , &std::log2  )
        ("erf"   , &std::erf   )
        ("erfc"  , &std::erfc  )
        ("lgamma", &std::lgamma)
        ("tgamma", &std::tgamma)
        ("trunc" , &std::trunc )
        ("round" , &std::round )
        //      ("crbt"  , &std::crbt  )
        ;

      binary_function.add
        ("pow"  , &std::pow      )
        ("min"  , &parser::my_min)
        ("max"  , &parser::my_max)
        ("atan2", &std::atan2    )
        ("fmod" , &std::fmod     )
        ("hypot", &std::hypot    )
          ;

#if !defined(AKANTU_NDEBUG)
      phx::function<algebraic_error_handler_> const error_handler = algebraic_error_handler_();
      qi::on_error<qi::fail>(start, error_handler(lbs::_4, lbs::_3, lbs::_2));
#endif

      expr    .name("expression");
      term    .name("term");
      factor  .name("factor");
      number  .name("numerical-value");
      variable.name("variable");
      function.name("function");

      constant.name("constants-list");
      unary_function.name("unary-functions-list");
      binary_function.name("binary-functions-list");

#if !defined AKANTU_NDEBUG
      if(AKANTU_DEBUG_TEST(dblDebug)) {
        qi::debug(expr);
        qi::debug(term);
        qi::debug(factor);
        qi::debug(number);
        qi::debug(variable);
        qi::debug(function);
      }
#endif
    }
    /* clang-format on */
  private:
    qi::rule<Iterator, Real(), Skipper> start;
    qi::rule<Iterator, Real(), Skipper> expr;
    qi::rule<Iterator, Real(), Skipper> term;
    qi::rule<Iterator, Real(), Skipper> factor;
    qi::rule<Iterator, Real(), Skipper> number;
    qi::rule<Iterator, Real(), Skipper> variable;
    qi::rule<Iterator, Real(), Skipper> function;

    qi::rule<Iterator, std::string(), Skipper> key;

    qi::real_parser<Real, qi::real_policies<Real>> real;

    qi::symbols<char, Real> constant;
    qi::symbols<char, Real (*)(Real)> unary_function;
    qi::symbols<char, Real (*)(Real, Real)> binary_function;

    const ParserSection & section;
  };

  /* ---------------------------------------------------------------------- */
  /* Vector Parser                                                          */
  /* ---------------------------------------------------------------------- */
  struct parsable_vector {
    operator Vector<Real>() {
      Vector<Real> tmp(_cells.size());
      auto it = _cells.begin();
      for (UInt i = 0; it != _cells.end(); ++it, ++i)
        tmp(i) = *it;
      return tmp;
    }

    std::vector<Real> _cells;
  };

  inline std::ostream & operator<<(std::ostream & stream,
                                   const parsable_vector & pv) {
    stream << "pv[";
    auto it = pv._cells.begin();
    if (it != pv._cells.end()) {
      stream << *it;
      for (++it; it != pv._cells.end(); ++it)
        stream << ", " << *it;
    }
    stream << "]";
    return stream;
  }

  struct parsable_matrix {
    operator Matrix<Real>() {
      size_t cols = 0;
      auto it_rows = _cells.begin();
      for (; it_rows != _cells.end(); ++it_rows)
        cols = std::max(cols, it_rows->_cells.size());

      Matrix<Real> tmp(_cells.size(), _cells[0]._cells.size(), 0.);

      it_rows = _cells.begin();
      for (UInt i = 0; it_rows != _cells.end(); ++it_rows, ++i) {
        auto it_cols = it_rows->_cells.begin();
        for (UInt j = 0; it_cols != it_rows->_cells.end(); ++it_cols, ++j) {
          tmp(i, j) = *it_cols;
        }
      }
      return tmp;
    }

    std::vector<parsable_vector> _cells;
  };

  inline std::ostream & operator<<(std::ostream & stream,
                                   const parsable_matrix & pm) {
    stream << "pm[";
    auto it = pm._cells.begin();
    if (it != pm._cells.end()) {
      stream << *it;
      for (++it; it != pm._cells.end(); ++it)
        stream << ", " << *it;
    }
    stream << "]";
    return stream;
  }

  /* ---------------------------------------------------------------------- */
  template <typename T1, typename T2>
  static void cont_add(T1 & cont, T2 & value) {
    cont._cells.push_back(value);
  }

  /* ---------------------------------------------------------------------- */
  template <class Iterator, typename Skipper = spirit::unused_type>
  struct VectorGrammar : qi::grammar<Iterator, parsable_vector(), Skipper> {
    VectorGrammar(const ParserSection & section)
        : VectorGrammar::base_type(start, "vector_algebraic_grammar"),
          number(section) {

      start = '[' > vector > ']';

      vector =
          (number[phx::bind(&cont_add<parsable_vector, Real>, lbs::_a,
                            lbs::_1)] >>
           *(',' >> number[phx::bind(&cont_add<parsable_vector, Real>, lbs::_a,
                                     lbs::_1)]))[lbs::_val = lbs::_a];

#if !defined(AKANTU_NDEBUG)
      phx::function<algebraic_error_handler_> const error_handler =
          algebraic_error_handler_();
      qi::on_error<qi::fail>(start, error_handler(lbs::_4, lbs::_3, lbs::_2));
#endif

      start.name("start");
      vector.name("vector");
      number.name("value");

#if !defined AKANTU_NDEBUG
      if (AKANTU_DEBUG_TEST(dblDebug)) {
        qi::debug(start);
        qi::debug(vector);
      }
#endif
    }

  private:
    qi::rule<Iterator, parsable_vector(), Skipper> start;
    qi::rule<Iterator, parsable_vector(), qi::locals<parsable_vector>, Skipper>
        vector;
    qi::rule<Iterator, Real(), Skipper> value;
    AlgebraicGrammar<Iterator, Skipper> number;
  };

  /* ---------------------------------------------------------------------- */
  static inline bool vector_eval(const ID & a, const ParserSection & section,
                                 parsable_vector & result) {
    std::string value = section.getParameter(a, _ppsc_current_and_parent_scope);
    std::string::const_iterator b = value.begin();
    std::string::const_iterator e = value.end();
    parser::VectorGrammar<std::string::const_iterator, qi::space_type> grammar(
        section);
    return qi::phrase_parse(b, e, grammar, qi::space, result);
  }

  /* ---------------------------------------------------------------------- */
  template <class Iterator, typename Skipper = spirit::unused_type>
  struct MatrixGrammar : qi::grammar<Iterator, parsable_matrix(), Skipper> {
    MatrixGrammar(const ParserSection & section)
        : MatrixGrammar::base_type(start, "matrix_algebraic_grammar"),
          vector(section) {

      start = '[' >> matrix >> ']';

      matrix =
          (rows[phx::bind(&cont_add<parsable_matrix, parsable_vector>, lbs::_a,
                          lbs::_1)] >>
           *(',' >> rows[phx::bind(&cont_add<parsable_matrix, parsable_vector>,
                                   lbs::_a, lbs::_1)]))[lbs::_val = lbs::_a];

      rows = eval_vector | vector;

      eval_vector = (key[lbs::_pass = phx::bind(&vector_eval, lbs::_1, section,
                                                lbs::_a)])[lbs::_val = lbs::_a];

      key = qi::char_("a-zA-Z_") >>
            *qi::char_("a-zA-Z_0-9") // coming from the InputFileGrammar
          ;

#if !defined(AKANTU_NDEBUG)
      phx::function<algebraic_error_handler_> const error_handler =
          algebraic_error_handler_();
      qi::on_error<qi::fail>(start, error_handler(lbs::_4, lbs::_3, lbs::_2));
#endif

      start.name("matrix");
      matrix.name("all_rows");
      rows.name("rows");
      vector.name("vector");
      eval_vector.name("eval_vector");

#ifndef AKANTU_NDEBUG
      if (AKANTU_DEBUG_TEST(dblDebug)) {
        qi::debug(start);
        qi::debug(matrix);
        qi::debug(rows);
        qi::debug(eval_vector);
        qi::debug(key);
      }
#endif
    }

  private:
    qi::rule<Iterator, parsable_matrix(), Skipper> start;
    qi::rule<Iterator, parsable_matrix(), qi::locals<parsable_matrix>, Skipper>
        matrix;
    qi::rule<Iterator, parsable_vector(), Skipper> rows;
    qi::rule<Iterator, parsable_vector(), qi::locals<parsable_vector>, Skipper>
        eval_vector;
    qi::rule<Iterator, std::string(), Skipper> key;
    VectorGrammar<Iterator, Skipper> vector;
  };

  /* ---------------------------------------------------------------------- */
  /* Randon Generator                                                       */
  /* ---------------------------------------------------------------------- */
  struct ParsableRandomGenerator {
    ParsableRandomGenerator(
        Real base = Real(),
        const RandomDistributionType & type = _rdt_not_defined,
        const parsable_vector & parameters = parsable_vector())
        : base(base), type(type), parameters(parameters) {}

    Real base;
    RandomDistributionType type;
    parsable_vector parameters;
  };

  inline std::ostream & operator<<(std::ostream & stream,
                                   const ParsableRandomGenerator & prg) {
    stream << "prg[" << prg.base << " " << UInt(prg.type) << " "
           << prg.parameters << "]";
    return stream;
  }

  /* ---------------------------------------------------------------------- */
  template <class Iterator, typename Skipper = spirit::unused_type>
  struct RandomGeneratorGrammar
      : qi::grammar<Iterator, ParsableRandomGenerator(), Skipper> {
    RandomGeneratorGrammar(const ParserSection & section)
        : RandomGeneratorGrammar::base_type(start, "random_generator_grammar"),
          number(section) {

      start = generator.alias();

      generator =
          qi::hold[distribution[lbs::_val = lbs::_1]] |
          number[lbs::_val = phx::construct<ParsableRandomGenerator>(lbs::_1)];

      distribution = (number >> generator_type >> '[' >> generator_params >>
                      ']')[lbs::_val = phx::construct<ParsableRandomGenerator>(
                               lbs::_1, lbs::_2, lbs::_3)];

      generator_params =
          (number[phx::bind(&cont_add<parsable_vector, Real>, lbs::_a,
                            lbs::_1)] >>
           *(',' > number[phx::bind(&cont_add<parsable_vector, Real>, lbs::_a,
                                    lbs::_1)]))[lbs::_val = lbs::_a];

#define AKANTU_RANDOM_DISTRIBUTION_TYPE_ADD(r, data, elem)                     \
  (BOOST_PP_STRINGIZE(BOOST_PP_TUPLE_ELEM(2, 0, elem)),                        \
   AKANTU_RANDOM_DISTRIBUTION_TYPES_PREFIX(BOOST_PP_TUPLE_ELEM(2, 0, elem)))

      generator_type.add BOOST_PP_SEQ_FOR_EACH(
          AKANTU_RANDOM_DISTRIBUTION_TYPE_ADD, _,
          AKANTU_RANDOM_DISTRIBUTION_TYPES);
#undef AKANTU_RANDOM_DISTRIBUTION_TYPE_ADD

#if !defined(AKANTU_NDEBUG)
      phx::function<algebraic_error_handler_> const error_handler =
          algebraic_error_handler_();
      qi::on_error<qi::fail>(start, error_handler(lbs::_4, lbs::_3, lbs::_2));
#endif

      start.name("random-generator");
      generator.name("random-generator");
      distribution.name("random-distribution");
      generator_type.name("generator-type");
      generator_params.name("generator-parameters");
      number.name("number");

#ifndef AKANTU_NDEBUG
      if (AKANTU_DEBUG_TEST(dblDebug)) {
        qi::debug(generator);
        qi::debug(distribution);
        qi::debug(generator_params);
      }
#endif
    }

  private:
    qi::rule<Iterator, ParsableRandomGenerator(), Skipper> start;
    qi::rule<Iterator, ParsableRandomGenerator(), Skipper> generator;
    qi::rule<Iterator, ParsableRandomGenerator(), Skipper> distribution;
    qi::rule<Iterator, parsable_vector(), qi::locals<parsable_vector>, Skipper>
        generator_params;
    AlgebraicGrammar<Iterator, Skipper> number;
    qi::symbols<char, RandomDistributionType> generator_type;
  };
} // namespace parser
} // namespace akantu

#endif /* __AKANTU_ALGEBRAIC_PARSER_HH__ */
