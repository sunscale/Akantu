/**
 * @file   aka_named_argument.hh
 *
 * @author Marco Arena
 *
 * @date creation  Fri Jun 16 2017
 *
 * @brief A Documented file.
 *
 * @section LICENSE
 *
 * Public Domain ? https://gist.github.com/ilpropheta/7576dce4c3249df89f85
 *
 */
/* -------------------------------------------------------------------------- */
#include <tuple>
#include <type_traits>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_AKA_NAMED_ARGUMENT_HH__
#define __AKANTU_AKA_NAMED_ARGUMENT_HH__

namespace akantu {

struct use_named_args_t {};
extern use_named_args_t use_named_args;

namespace named_argument {
  /* -- Pack utils (proxy version) --------------------------------------------
   */
  /// Proxy containing [tag, value]
  template <typename tag, typename type> struct param_t {
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
      _tmp = (std::is_same<T, typename std::decay<head>::type::_tag>::value)
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

  template <typename T, typename head, typename... tail> struct type_at {
    enum { _tmp = type_at_p<T, head, tail...>::_pos };
    enum { _pos = _tmp == 1 ? 0 : (_tmp == -1 ? -1 : _tmp - 1) };
  };

  /*  Same as get_at but it's supposed to be used by passing
      a pack of param_t (_type is retrieved instead)
      This and get_at should be refactored.
  */
  template <int pos, typename T> struct get_required {
    static_assert(pos >= 0, " Missing required parameter, check the value of T "
                            "in get_required<pos, T>");

    template <typename... types>
    static decltype(auto) get(std::tuple<types &...> && t) {
      return std::forward<decltype(std::get<pos>(t)._value)>(
          std::get<pos>(t)._value);
    }
  };

  // Optional version
  template <int pos> struct get_optional {
    template <typename T, typename... types>
    static decltype(auto) get(T &&, std::tuple<types &...> && t) {
      return get_required<pos, T>::get(std::forward<decltype(t)>(t));
    }
  };

  template <> struct get_optional<-1> {
    template <typename T, typename... types>
    static decltype(auto) get(T && _default, std::tuple<types &...> &&) {
      return std::forward<T>(_default);
    }
  };

} // namespace named_argument

// CONVENIENCE MACROS FOR CLASS DESIGNERS ==========
#define TAG_OF_ARGUMENT(name) p_##name
#define TAG_OF_ARGUMENT_WNS(name) named_argument::TAG_OF_ARGUMENT(name)

#define REQUIRED_NAMED_ARG(name)                                               \
  named_argument::get_required<                                                \
      named_argument::type_at<TAG_OF_ARGUMENT_WNS(name), pack...>::_pos,       \
      TAG_OF_ARGUMENT_WNS(                                                     \
          name)>::get(std::forward_as_tuple<pack...>(_pack...))

#define OPTIONAL_NAMED_ARG(name, _defaultVal)                                  \
  named_argument::get_optional<                                                \
      named_argument::type_at<TAG_OF_ARGUMENT_WNS(name), pack...>::_pos>::     \
      get(_defaultVal, std::forward_as_tuple<pack...>(_pack...))

#define DECLARE_NAMED_ARGUMENT(name)                                           \
  namespace named_argument {                                                   \
    struct TAG_OF_ARGUMENT(name) {};                                           \
  }                                                                            \
  extern named_argument::param_proxy<TAG_OF_ARGUMENT_WNS(name)> _##name

#define CREATE_NAMED_ARGUMENT(name)                                            \
  named_argument::param_proxy<TAG_OF_ARGUMENT_WNS(name)> _##name

// FULL EXAMPLE ====================================
// step 1) generate tags you need
// namespace {
// CREATE_TAG(title);
// CREATE_TAG(h);
// CREATE_TAG(w);
// CREATE_TAG(posx);
// CREATE_TAG(posy);
// CREATE_TAG(handle);
// } // namespace

// // step 2) use tags/proxies in your class
// class window {
// public:
//   // classic constructor
//   window(string pTitle, int pH, int pW, int pPosx, int pPosy, int &
//   pHandle)
//       : title(move(pTitle)), h(pH), w(pW), posx(pPosx), posy(pPosy),
//         handle(pHandle) {}

//   // constructor using proxies (e.g. _title = "title")
//   template <typename... pack>
//   window(use_named_t, pack &&... _pack)
//       : window{REQUIRED_NAME(title),   // required
//                OPTIONAL_NAME(h, 100),  // optional
//                OPTIONAL_NAME(w, 400),  // optional
//                OPTIONAL_NAME(posx, 0), // optional
//                OPTIONAL_NAME(posy, 0), // optional
//                REQUIRED_NAME(handle)}  // required
//   {}

//   // constructor using tags (e.g. __title, "title")
//   template <typename... pack>
//   window(use_tags_t, pack &&... _pack)
//       : window{REQUIRED_TAG(title),   // required
//                OPTIONAL_TAG(h, 100),  // optional
//                OPTIONAL_TAG(w, 400),  // optional
//                OPTIONAL_TAG(posx, 0), // optional
//                OPTIONAL_TAG(posy, 0), // optional
//                REQUIRED_TAG(handle)}  // required
//   {}

//   friend ostream & operator<<(ostream & os, const window &);

// private:
//   string title;
//   int h, w;
//   int posx, posy;
//   int & handle;
// };

// ostream & operator<<(ostream & os, const window & w) {
//   os << "WINDOW: " << w.title << " (" << w.h << "x" << w.w << "), at " <<
//   w.posx
//      << "," << w.posy << "->" << w.handle;
//   return os;
// }

// int main() {
//   int i = 5;
//   {
//     window w{use_tags, __title, "Title", __h, 10, __w, 100, __handle, i};
//     cout << w << endl;
//   }

//   {
//     window w{use_named, _h = 10, _title = "Title", _handle = i, _w = 100};
//     cout << w << endl;
//   }

//   {
//     window w{"Title", 10, 400, 0, 0, i};
//     cout << w << endl;
//   }
// }

} // namespace akantu

#endif /* __AKANTU_AKA_NAMED_ARGUMENT_HH__ */
