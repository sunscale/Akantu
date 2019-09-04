/**
 * @file   input_file_parser.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Nov 13 2013
 * @date last modification: Mon Dec 18 2017
 *
 * @brief  Grammar definition for the input files
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
// Boost
/* -------------------------------------------------------------------------- */

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/phoenix_bind.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/variant/recursive_variant.hpp>

#ifndef __AKANTU_INPUT_FILE_PARSER_HH__
#define __AKANTU_INPUT_FILE_PARSER_HH__

namespace spirit = boost::spirit;
namespace qi = boost::spirit::qi;
namespace lbs = boost::spirit::qi::labels;
namespace ascii = boost::spirit::ascii;
namespace phx = boost::phoenix;

namespace akantu {

namespace parser {
  struct error_handler_ {
    template <typename, typename, typename, typename> struct result {
      using type = void;
    };

    template <typename Iterator>
    void operator()(qi::info const & what, Iterator err_pos, Iterator /*first*/,
                    Iterator /*last*/) const {
      spirit::classic::file_position pos = err_pos.get_position();

      AKANTU_EXCEPTION("Parse error [ "
                       << "Expecting " << what << " instead of \"" << *err_pos
                       << "\" ]"
                       << " in file " << pos.file << " line " << pos.line
                       << " column " << pos.column << std::endl
                       << "'" << err_pos.get_currentline() << "'" << std::endl
                       << std::setw(pos.column) << " "
                       << "^- here");
    }

  private:
  };

  static ParserSection & create_subsection(
      const ParserType & type, const boost::optional<std::string> & opt_name,
      const boost::optional<std::string> & opt_option, ParserSection & sect) {
    std::string option = "";
    if (opt_option)
      option = *opt_option;

    static size_t id = 12;
    std::string name = "anonymous_" + std::to_string(id++);
    if (opt_name)
      name = *opt_name;

    ParserSection sect_tmp(name, type, option, sect);
    return sect.addSubSection(sect_tmp);
  }

  template <typename Iter>
  static bool create_parameter(boost::iterator_range<Iter> & rng,
                               std::string & value, ParserSection & sect) {
    try {
      std::string name(rng.begin(), rng.end());
      name = trim(name);
      spirit::classic::file_position pos = rng.begin().get_position();

      ParserParameter param_tmp(name, value, sect);
      param_tmp.setDebugInfo(pos.file, pos.line, pos.column);
      sect.addParameter(param_tmp);
    } catch (debug::Exception & e) {
      return false;
    }
    return true;
  }

  static std::string concatenate(const std::string & t1,
                                 const std::string & t2) {
    return (t1 + t2);
  }

  /* ---------------------------------------------------------------------- */
  /* Grammars definitions                                                   */
  /* ---------------------------------------------------------------------- */
  template <class Iterator>
  struct InputFileGrammar
      : qi::grammar<Iterator, void(), typename Skipper<Iterator>::type> {
    InputFileGrammar(ParserSection * sect)
        : InputFileGrammar::base_type(start, "input_file_grammar"),
          parent_section(sect) {

      /* clang-format off */
    start
      =   mini_section(parent_section)
      ;

    mini_section
      =   *(
               entry  (lbs::_r1)
           |   section(lbs::_r1)
           )
      ;

    entry
      =   (
             qi::raw[key]
          >> '='
          >  value
          ) [ lbs::_pass = phx::bind(&create_parameter<Iterator>,
                                     lbs::_1,
                                     lbs::_2,
                                     *lbs::_r1) ]
      ;

    section
      =   (
            qi::no_case[section_type]
          > qi::lexeme
               [
                   -section_name
                   > -section_option
               ]
           ) [ lbs::_a = &phx::bind(&create_subsection,
                                    lbs::_1,
                                    phx::at_c<0>(lbs::_2),
                                    phx::at_c<1>(lbs::_2),
                                    *lbs::_r1) ]
          > '['
          > mini_section(lbs::_a)
          > ']'
      ;

    section_name
      =   qi::char_("a-zA-Z_") >> *qi::char_("a-zA-Z_0-9")
      ;

    section_option
      =   (+ascii::space >> section_name) [ lbs::_val = lbs::_2 ]
      ;

    key
      =   qi::char_("a-zA-Z_") >> *qi::char_("a-zA-Z_0-9")
      ;

    value
      =   (
           mono_line_value          [ lbs::_a = phx::bind(&concatenate, lbs::_a, lbs::_1) ]
        > *(
            '\\' > mono_line_value  [ lbs::_a = phx::bind(&concatenate, lbs::_a, lbs::_1) ]
           )
          ) [ lbs::_val = lbs::_a ]
      ;

    mono_line_value
      =   qi::lexeme
             [
           +(qi::char_ - (qi::char_('=') | spirit::eol | '#' | ';' | '\\'))
             ]
      ;

    skipper
      =     ascii::space
      |     "#" >> *(qi::char_ - spirit::eol)
      ;

      /* clang-format on */

#define AKANTU_SECTION_TYPE_ADD(r, data, elem)                                 \
  (BOOST_PP_STRINGIZE(elem), BOOST_PP_CAT(ParserType::_, elem))

      section_type.add BOOST_PP_SEQ_FOR_EACH(AKANTU_SECTION_TYPE_ADD, _,
                                             AKANTU_SECTION_TYPES);
#undef AKANTU_SECTION_TYPE_ADD

#if !defined(AKANTU_NDEBUG)
      phx::function<error_handler_> const error_handler = error_handler_();
      qi::on_error<qi::fail>(start,
                             error_handler(lbs::_4, lbs::_3, lbs::_1, lbs::_2));
#endif

      section.name("section");
      section_name.name("section-name");
      section_option.name("section-option");
      mini_section.name("section-content");
      entry.name("parameter");
      key.name("parameter-name");
      value.name("parameter-value");
      section_type.name("section-types-list");
      mono_line_value.name("mono-line-value");

#if !defined AKANTU_NDEBUG
      if (AKANTU_DEBUG_TEST(dblDebug)) {
        //        qi::debug(section);
        qi::debug(section_name);
        qi::debug(section_option);
        //        qi::debug(mini_section);
        //	qi::debug(entry);
        qi::debug(key);
        qi::debug(value);
        qi::debug(mono_line_value);
      }
#endif
    }

    const std::string & getErrorMessage() const { return error_message; };

    using skipper_type = typename Skipper<Iterator>::type;
    skipper_type skipper;

  private:
    std::string error_message;

    qi::rule<Iterator, void(ParserSection *), skipper_type> mini_section;
    qi::rule<Iterator, void(ParserSection *), qi::locals<ParserSection *>,
             skipper_type>
        section;
    qi::rule<Iterator, void(), skipper_type> start;
    qi::rule<Iterator, std::string()> section_name;
    qi::rule<Iterator, std::string()> section_option;
    qi::rule<Iterator, void(ParserSection *), skipper_type> entry;
    qi::rule<Iterator, std::string(), skipper_type> key;
    qi::rule<Iterator, std::string(), qi::locals<std::string>, skipper_type>
        value;
    qi::rule<Iterator, std::string(), skipper_type> mono_line_value;

    qi::symbols<char, ParserType> section_type;

    ParserSection * parent_section;
  };
} // namespace parser

} // namespace akantu

#endif /* __AKANTU_INPUT_FILE_PARSER_HH__ */
