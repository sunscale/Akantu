/* -------------------------------------------------------------------------- */
#include "aka_array.hh"
/* -------------------------------------------------------------------------- */
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
/* -------------------------------------------------------------------------- */

namespace py = pybind11;
namespace _aka = akantu;

namespace akantu {

template <typename VecType, typename T> class Proxy : public VecType {

public:
  Proxy(T * data, UInt size) {
    this->values = data;
    this->wrapped = true;
    this->n[0] = size;
  }

  Proxy(T * data, UInt size1, UInt size2) {
    this->values = data;
    this->wrapped = true;
    this->n[0] = size1;
    this->n[1] = size2;
  }

  ~Proxy() { this->values = nullptr; }
};

template <typename T> class Proxy<_aka::Array<T>, T> : public _aka::Array<T> {
protected:
  // deallocate the memory
  void deallocate() override final {}

  // allocate the memory
  void allocate(__attribute__((unused)) UInt size,
                __attribute__((unused)) UInt nb_component) override final {}

  // allocate and initialize the memory
  void allocate(__attribute__((unused)) UInt size,
                __attribute__((unused)) UInt nb_component,
                __attribute__((unused)) const T & value) override final {}

public:
  Proxy(T * data, UInt size, UInt nb_component) {
    this->values = data;
    this->size_ = size;
    this->nb_component = nb_component;
  }

  Proxy(const Array<T> & src) {
    this->values = src.storage();
    this->size_ = src.size();
    this->nb_component = src.getNbComponent();
  }

  ~Proxy() { this->values = nullptr; }

  void resize(__attribute__((unused)) UInt size,
              __attribute__((unused)) const T & val) override final {
    AKANTU_EXCEPTION("cannot resize a temporary array");
  }

  void resize(__attribute__((unused)) UInt new_size) override final {
    AKANTU_EXCEPTION("cannot resize a temporary array");
  }

  void reserve(__attribute__((unused)) UInt new_size) override final {
    AKANTU_EXCEPTION("cannot resize a temporary array");
  }
};

template <typename T> using vec_proxy = Proxy<Vector<T>, T>;
template <typename T> using mat_proxy = Proxy<Matrix<T>, T>;
template <typename T> using array_proxy = Proxy<Array<T>, T>;

} // namespace akantu

namespace pybind11 {
namespace detail {

  template <typename U>
  using array_type = array_t<U, array::c_style | array::forcecast>;

  template <typename T>
  void create_proxy(std::unique_ptr<_aka::vec_proxy<T>> & proxy,
                    array_type<T> ref) {
    proxy =
        std::make_unique<_aka::vec_proxy<T>>(ref.mutable_data(), ref.shape(0));
  }

  template <typename T>
  void create_proxy(std::unique_ptr<_aka::mat_proxy<T>> & proxy,
                    array_type<T> ref) {
    proxy = std::make_unique<_aka::mat_proxy<T>>(ref.mutable_data(),
                                                 ref.shape(0), ref.shape(1));
  }

  template <typename T>
  void create_proxy(std::unique_ptr<_aka::array_proxy<T>> & proxy,
                    array_type<T> ref) {
    proxy = std::make_unique<_aka::array_proxy<T>>(ref.mutable_data(),
                                                   ref.shape(0), ref.shape(1));
  }

  /* ------------------------------------------------------------------------ */
  template <typename T>
  py::handle aka_array_cast(const _aka::Array<T> & src,
                            py::handle base = handle(), bool writeable = true) {
    array a;
    a = array_type<T>({src.size(), src.getNbComponent()}, src.storage(), base);

    if (not writeable)
      array_proxy(a.ptr())->flags &= ~detail::npy_api::NPY_ARRAY_WRITEABLE_;

    return a.release();
  }

  template <typename T>
  py::handle aka_array_cast(const _aka::Vector<T> & src,
                            py::handle base = handle(), bool writeable = true) {
    array a;
    a = array_type<T>({src.size()}, src.storage(), base);

    if (not writeable)
      array_proxy(a.ptr())->flags &= ~detail::npy_api::NPY_ARRAY_WRITEABLE_;

    return a.release();
  }

  template <typename T>
  py::handle aka_array_cast(const _aka::Matrix<T> & src,
                            py::handle base = handle(), bool writeable = true) {
    array a;
    a = array_type<T>({src.size(0), src.size(1)}, src.storage(), base);

    if (not writeable)
      array_proxy(a.ptr())->flags &= ~detail::npy_api::NPY_ARRAY_WRITEABLE_;

    return a.release();
  }

  /* ------------------------------------------------------------------------ */
  template <typename VecType, typename T>
  class [[gnu::visibility("default")]] my_type_caster  {
  protected:
    using type = VecType;
    using proxy_type = typename _aka::Proxy<VecType, T>;
    type value;

  public:
    static PYBIND11_DESCR name() { return type_descr(_("Toto")); };

    /**
     * Conversion part 1 (Python->C++)
     */
    bool load(handle src, bool convert) {
      bool need_copy = not isinstance<array_type<T>>(src);

      auto && fits = [&](auto && aref) {
        auto && dims = aref.ndim();
        if (dims < 1 || dims > 2)
          return false;

        return true;
      };

      if (not need_copy) {
        // We don't need a converting copy, but we also need to check whether
        // the strides are compatible with the Ref's stride requirements
        auto aref = py::cast<array_type<T>>(src);

        if (not fits(aref)) {
          return false;
        }
        copy_or_ref = std::move(aref);
      } else {
        if (not convert) {
          return false;
        }

        auto copy = array_type<T>::ensure(src);
        if (not copy) {
          return false;
        }

        if (not fits(copy)) {
          return false;
        }
        copy_or_ref = std::move(array_type<T>::ensure(src));
        loader_life_support::add_patient(copy_or_ref);
      }

      create_proxy(array_proxy, copy_or_ref);
      return true;
    }

    operator type *() { return array_proxy.get(); }
    operator type &() { return *array_proxy; }

    template <typename _T>
    using cast_op_type = pybind11::detail::cast_op_type<_T>;

    /**
     * Conversion part 2 (C++ -> Python)
     */
    static handle cast(const type & src, return_value_policy policy,
                       handle parent) {
      switch (policy) {
      case return_value_policy::copy:
        return aka_array_cast<T>(src);
      case return_value_policy::reference_internal:
        return aka_array_cast<T>(src, parent);
      case return_value_policy::reference:
      case return_value_policy::automatic:
      case return_value_policy::automatic_reference:
        return aka_array_cast<T>(src, none());
      default:
        pybind11_fail("Invalid return_value_policy for ArrayProxy type");
      }
    }

  protected:
    std::unique_ptr<proxy_type> array_proxy;
    array_type<T> copy_or_ref;
  };

  /* ------------------------------------------------------------------------ */
  // specializations
  /* ------------------------------------------------------------------------ */

  template <typename T>
  struct type_caster<_aka::Array<T>>
      : public my_type_caster<_aka::Array<T>, T> {
    using my_type_caster<_aka::Array<T>, T>::cast_op_type;
    using my_type_caster<_aka::Array<T>, T>::cast;
    using my_type_caster<_aka::Array<T>, T>::load;
  };

  template <typename T>
  struct type_caster<_aka::Vector<T>>
      : public my_type_caster<_aka::Vector<T>, T> {
    using my_type_caster<_aka::Vector<T>, T>::cast_op_type;
    using my_type_caster<_aka::Vector<T>, T>::cast;
    using my_type_caster<_aka::Vector<T>, T>::load;
  };

  template <typename T>
  struct type_caster<_aka::Matrix<T>>
      : public my_type_caster<_aka::Matrix<T>, T> {
    using my_type_caster<_aka::Matrix<T>, T>::cast_op_type;
    using my_type_caster<_aka::Matrix<T>, T>::cast;
    using my_type_caster<_aka::Matrix<T>, T>::load;
  };

  // template <typename T>
  // struct type_caster<_aka::Vector<T>>
  //     : public my_type_caster<_aka::Vector<T>, T> {};

  // /**
  //  * Type caster for Vector classes
  //  */
  // template <template <typename> class V, typename T> struct type_caster<V<T>>
  // {
  //   using type = V<T>;

  // public:
  //   PYBIND11_TYPE_CASTER(type, _("Vector<T>"));

  //   /**
  //    * Conversion part 1 (Python->C++): convert a PyObject into a Vector
  //    * instance or return false upon failure. The second argument
  //    * indicates whether implicit conversions should be applied.
  //    */
  //   bool load(handle src, bool convert) {
  //     if (!convert && !array_type<typename type::value_type>::check_(src))
  //       return false;

  //     auto buf = array_type<typename type::value_type>::ensure(src);
  //     value.move(_aka::VectorProxy<T>(buf));

  //     return true;
  //   }

  //   /**
  //    * Conversion part 2 (C++ -> Python): convert a grid instance into
  //    * a Python object. The second and third arguments are used to
  //    * indicate the return value policy and parent object (for
  //    * ``return_value_policy::reference_internal``) and are generally
  //    * ignored by implicit casters.
  //    *
  //    * TODO: do not ignore policy (see pybind11/eigen.h)
  //    */
  //   static handle cast(const type & src, return_value_policy /* policy */,
  //                      handle /*parent*/) {
  //     // none() passed as parent to get a correct nocopy
  //     auto a = array_type<typename type::value_type>(src.size(),
  //     src.storage(),
  //                                                    none());
  //     return a.release();
  //   }
  // };

} // namespace detail
} // namespace pybind11
