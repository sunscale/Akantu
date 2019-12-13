// Copyright (c) 2016 Vittorio Romeo
// License: AFL 3.0 | https://opensource.org/licenses/AFL-3.0
// http://vittorioromeo.info | vittorio.romeo@outlook.com

#ifndef AKANTU_AKA_STATIC_IF_HH
#define AKANTU_AKA_STATIC_IF_HH

#ifndef AKANTU_ITERATORS_NAMESPACE
#define AKANTU_ITERATORS_NAMESPACE akantu
#endif

#include <utility>

#define FWD(...) ::std::forward<decltype(__VA_ARGS__)>(__VA_ARGS__)

namespace AKANTU_ITERATORS_NAMESPACE {

template <typename TPredicate> auto static_if(TPredicate /*unused*/) noexcept;

namespace impl {
  template <bool TPredicateResult> struct static_if_impl;

  template <typename TFunctionToCall> struct static_if_result;

  template <typename TF> auto make_static_if_result(TF && f) noexcept;

  template <> struct static_if_impl<true> {
    template <typename TF> auto & else_(TF && /*unused*/) noexcept {
      // Ignore `else_`, as the predicate is true.
      return *this;
    }

    template <typename TPredicate>
    auto & else_if(TPredicate /*unused*/) noexcept {
      // Ignore `else_if`, as the predicate is true.
      return *this;
    }

    template <typename TF> auto then(TF && f) noexcept {
      // We found a matching branch, just make a result and
      // ignore everything else.
      return make_static_if_result(FWD(f));
    }
  };

  template <> struct static_if_impl<false> {
    template <typename TF> auto & then(TF && /*unused*/) noexcept {
      // Ignore `then`, as the predicate is false.
      return *this;
    }

    template <typename TF> auto else_(TF && f) noexcept {
      // (Assuming that `else_` is after all `else_if` calls.)

      // We found a matching branch, just make a result and
      // ignore everything else.

      return make_static_if_result(FWD(f));
    }

    template <typename TPredicate>
    auto else_if(TPredicate /*unused*/) noexcept {
      return static_if(TPredicate{});
    }

    template <typename... Ts> auto operator()(Ts &&... /*unused*/) noexcept {
      // If there are no `else` branches, we must ignore calls
      // to a failed `static_if` matching.
    }
  };

  template <typename TFunctionToCall>
  struct static_if_result : TFunctionToCall {
    // Perfect-forward the function in the result instance.
    template <typename TFFwd>
    explicit static_if_result(TFFwd && f) noexcept // NOLINT
        : TFunctionToCall(FWD(f)) {}

    // Ignore everything, we found a result.
    template <typename TF> auto & then(TF && /*unused*/) noexcept {
      return *this;
    }

    template <typename TPredicate>
    auto & else_if(TPredicate /*unused*/) noexcept {
      return *this;
    }

    template <typename TF> auto & else_(TF && /*unused*/) noexcept {
      return *this;
    }
  };

  template <typename TF> auto make_static_if_result(TF && f) noexcept {
    return static_if_result<TF>{FWD(f)};
  }
} // namespace impl

template <typename TPredicate> auto static_if(TPredicate /*unused*/) noexcept {
  return impl::static_if_impl<TPredicate{}>{};
}

#undef FWD
} // namespace AKANTU_ITERATORS_NAMESPACE

#endif /* AKANTU_AKA_STATIC_IF_HH */
