/* -------------------------------------------------------------------------- */
#include "aka_array.hh"
/* -------------------------------------------------------------------------- */
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
/* -------------------------------------------------------------------------- */

namespace py = pybind11;
namespace _aka = akantu;

namespace akantu {
template <typename T> class ArrayProxy : public Array<T> {
protected:
  // deallocate the memory
  void deallocate() override final {}

  // allocate the memory
  void allocate(UInt size, UInt nb_component) override final {}

  // allocate and initialize the memory
  void allocate(UInt size, UInt nb_component, const T & value) override final {}

public:
  ArrayProxy(T * data, UInt size, UInt nb_component) {
    this->values = data;
    this->size_ = size;
    this->nb_component = nb_component;
  }

  ArrayProxy(const Array<T> & src) {
    this->values = src.storage();
    this->size_ = src.size();
    this->nb_component = src.getNbComponent();
  }

  void resize(UInt size, const T & val) override final {
    AKANTU_EXCEPTION("cannot resize a temporary array");
  }

  void resize(UInt new_size) override final {
    AKANTU_EXCEPTION("cannot resize a temporary array");
  }

  void reserve(UInt new_size) override final {
    AKANTU_EXCEPTION("cannot resize a temporary array");
  }
};

} // namespace akantu

namespace pybind11 {
namespace detail {
  /* ------------------------------------------------------------------------ */
  template <typename T>
  py::handle aka_array_cast(const _aka::Array<T> & src,
                            py::handle base = handle(), bool writeable = true) {
    constexpr ssize_t elem_size = sizeof(T);
    array a;
    a = array({src.size(), src.getNbComponent()},
              {elem_size * src.getNbComponent(), elem_size}, src.storage(),
              base);

    if (not writeable)
      array_proxy(a.ptr())->flags &= ~detail::npy_api::NPY_ARRAY_WRITEABLE_;

    return a.release();
  }

  template <typename T> struct type_caster<_aka::Array<T>> {
  protected:
    using py_array = array_t<T, array::forcecast | array::c_style>;

  public:
    static constexpr auto name =
        _("numpy.ndarray[") + npy_format_descriptor<T>::name +
        _(", flags.writeable") + _(", flags.c_contiguous") + _("]");

    /**
     * Conversion part 1 (Python->C++)
     */
    bool load(handle src, bool convert) {
      bool need_copy = not isinstance<py_array>(src);

      auto && fits = [&](auto && a) {
        auto && dims = copy_or_ref.ndim();
        if (dims < 1 || dims > 2)
          return false;

        return true;
      };

      if (not need_copy) {
        // We don't need a converting copy, but we also need to check whether
        // the strides are compatible with the Ref's stride requirements
        py_array aref = reinterpret_borrow<py_array>(src);

        if (not fits(aref)) {
          return false;
        }

        copy_or_ref = std::move(aref);
      } else {
        if (not convert) {
          return false;
        }

        auto copy = py_array::ensure(src);
        if (not copy) {
          return false;
        }

        if (not fits(copy)) {
          return false;
        }
        copy_or_ref = std::move(py_array::ensure(src));
        loader_life_support::add_patient(copy_or_ref);
      }

      array_proxy.reset(new _aka::ArrayProxy<T>> (copy_or_ref.mutable_data(),
                                                  copy_or_ref.shape(0),
                                                  copy_or_ref.shape(1)));
      return true;
    }

    operator _aka::Array<T> *() { return array_proxy.get(); }
    operator _aka::Array<T> &() { return *array_proxy; }
    template <typename _T>
    using cast_op_type = pybind11::detail::cast_op_type<_T>;

    /**
     * Conversion part 2 (C++ -> Python)
     */
    static handle cast(const _aka::Array<T> & src, return_value_policy policy,
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
    std::unique_ptr<_aka::Array<T>> array_proxy;
    py_array copy_or_ref;
  };
} // namespace detail
} // namespace pybind11

namespace{
  py::module & register_arrays(py::module & mod) {

    py::class_<_aka::Vector<bool>>(mod, "VectorBool");
    py::class_<_aka::Vector<double>>(mod, "VectorReal");    
    return mod;
  }
}
