/**
 * @file   aka_named_argument.hh
 *
 * @author Marco Arena
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 16 2017
 * @date last modification: Wed Dec 06 2017
 *
 * @brief  tool to use named arguments in functions
 *
 *
 * Public Domain ? https://gist.github.com/ilpropheta/7576dce4c3249df89f85
 *
 */
/* -------------------------------------------------------------------------- */
#include "aka_compatibilty_with_cpp_standard.hh"
/* -------------------------------------------------------------------------- */
#include <tuple>
#include <type_traits>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_AKA_NAMED_ARGUMENT_HH__
#define __AKANTU_AKA_NAMED_ARGUMENT_HH__

namespace akantu {

namespace named_argument {
  struct param_t_trait {};

  /* -- Pack utils (proxy version) ------------------------------------------ */
  /// Proxy containing [tag, value]
  template <typename tag, typename type> struct param_t : param_t_trait {
    using _tag = tag;
    using _type = type;

    template <typename T>
    explicit param_t(T && value) : _value(std::forward<T>(value)) {}

    type _value;
  };

  /*
   * Tagged proxy that allows syntax _name = value
   * operator=(T&&) returns a param_t instance
   **/
  template <typename tag> struct param_proxy {
    using _tag = tag;

    template <typename T> decltype(auto) operator=(T && value) {
      return param_t<tag, decltype(value)>{std::forward<T>(value)};
    }
  };

  /*  Same as type_at but it's supposed to be used by passing
      a pack of param_t (_tag is looked for instead of a
      plain type). This and type_at should be refactored.
  */
  template <typename T, typename head, typename... tail> struct type_at_p {
    enum {
      _tmp = (std::is_same<T, typename std::decay_t<head>::_tag>::value)
                 ? 0
                 : type_at_p<T, tail...>::_pos
    };
    enum { _pos = _tmp == -1 ? -1 : 1 + _tmp };
  };

  template <typename T, typename head> struct type_at_p<T, head> {
    enum {
      _pos =
          (std::is_same<T, typename std::decay<head>::type::_tag>::value ? 1
                                                                         : -1)
    };
  };

  template <typename... Ts> struct type_at {
    enum { _pos = -1 };
  };

  template <typename T, typename head, typename... tail>
  struct type_at<T, head, tail...> {
    enum { _tmp = type_at_p<T, head, tail...>::_pos };
    enum { _pos = _tmp == 1 ? 0 : (_tmp == -1 ? -1 : _tmp - 1) };
  };

  /*  Same as get_at but it's supposed to be used by passing
      a pack of param_t (_type is retrieved instead)
      This and get_at should be refactored.
  */
  template <int pos, int curr> struct get_at {
    static_assert(pos >= 0, "Required parameter");

    template <typename head, typename... tail>
    static decltype(auto) get(head &&, tail &&... t) {
      return get_at<pos, curr + 1>::get(std::forward<tail>(t)...);
    }
  };

  template <int pos> struct get_at<pos, pos> {
    static_assert(pos >= 0, "Required parameter");

    template <typename head, typename... tail>
    static decltype(auto) get(head && h, tail &&...) {
      return std::forward<decltype(h._value)>(h._value);
    }
  };

  // Optional version
  template <int pos, int curr> struct get_optional {
    template <typename T, typename... pack>
    static decltype(auto) get(T &&, pack &&... _pack) {
      return get_at<pos, curr>::get(std::forward<pack>(_pack)...);
    }
  };

  template <int curr> struct get_optional<-1, curr> {
    template <typename T, typename... pack>
    static decltype(auto) get(T && _default, pack &&...) {
      return std::forward<T>(_default);
    }
  };

} // namespace named_argument

// CONVENIENCE MACROS FOR CLASS DESIGNERS ==========
#define TAG_OF_ARGUMENT(_name) p_##_name
#define TAG_OF_ARGUMENT_WNS(_name) TAG_OF_ARGUMENT(_name)

#define REQUIRED_NAMED_ARG(_name)                                              \
  named_argument::get_at<                                                      \
      named_argument::type_at<TAG_OF_ARGUMENT_WNS(_name), pack...>::_pos,      \
      0>::get(std::forward<pack>(_pack)...)

#define REQUIRED_NAMED_ARG(_name)                                              \
  named_argument::get_at<                                                      \
      named_argument::type_at<TAG_OF_ARGUMENT_WNS(_name), pack...>::_pos,      \
      0>::get(std::forward<pack>(_pack)...)
#define OPTIONAL_NAMED_ARG(_name, _defaultVal)                                 \
  named_argument::get_optional<                                                \
      named_argument::type_at<TAG_OF_ARGUMENT_WNS(_name), pack...>::_pos,      \
      0>::get(_defaultVal, std::forward<pack>(_pack)...)

#define DECLARE_NAMED_ARGUMENT(name)                                           \
  struct TAG_OF_ARGUMENT(name) {};                                             \
  named_argument::param_proxy<TAG_OF_ARGUMENT_WNS(name)> _##name               \
      __attribute__((unused))

namespace {
  struct use_named_args_t {};
  use_named_args_t use_named_args __attribute__((unused));
} // namespace

template <typename T> struct is_named_argument : public std::false_type {};

template <typename... type>
struct is_named_argument<named_argument::param_t<type...>>
    : public std::true_type {};

template <typename... pack>
using are_named_argument =
    aka::conjunction<is_named_argument<std::decay_t<pack>>...>;

} // namespace akantu

#endif /* __AKANTU_AKA_NAMED_ARGUMENT_HH__ */
