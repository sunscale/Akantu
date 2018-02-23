#include "aka_array.hh"
#include "aka_types.hh"
/* -------------------------------------------------------------------------- */
#include <pybind11/cast.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_PYBIND11_AKANTU_HH__
#define __AKANTU_PYBIND11_AKANTU_HH__

namespace akantu {

template <typename T> class ArrayProxy : public Array<T> {
public:
  ArrayProxy(T * wrapped_memory, UInt size = 0, UInt nb_component = 1,
             const ID & id = "")
      : Array<T>(0, nb_component, id) {
    this->values = wrapped_memory;
    this->size_ = size;
  }

  ArrayProxy(const ArrayProxy<T> & array)
      : Array<T>(0, array.nb_component, array.id) {
    this->values = array.values;
    this->size_ = array.size_;
  }

  ArrayProxy(ArrayProxy<T> && array) {
    this->nb_component = std::move(array.nb_component);
    this->values = std::move(array.values);
    this->size_ = std::move(array.size_);
    this->id = std::move(array.id);
  }

  ArrayProxy(Array<T> & array)
      : Array<T>(0, array.getNbComponent(), array.getID()) {
    this->values = array.storage();
    this->size_ = array.size();
  }

  ~ArrayProxy() {
    this->values = nullptr;
    this->size_ = 0;
  }

  void setNbComponent(UInt nb_component) {
    UInt new_size = this->size_ / nb_component;
    AKANTU_DEBUG_ASSERT(
        nb_component * new_size == this->nb_component * this->size_,
        nb_component
            << " is not valid as a new number of component for this array");

    this->nb_component = nb_component;
    this->size_ = new_size;
  }

  void resize(UInt new_size) {
    AKANTU_DEBUG_ASSERT(this->size_ == new_size,
                        "cannot resize a temporary vector");
  }
};

template <typename T> decltype(auto) make_proxy(Array<T> & array) {
  return ArrayProxy<T>(array);
}

template <typename T> decltype(auto) make_proxy(const Matrix<T> & array) {
  return MatrixProxy<T>(array);
}

} // namespace akantu

// namespace pybind11 {
// namespace detail {
//   template <> struct type_caster<akantu::ArrayProxy<akantu::Real>> {
//   public:
//     PYBIND11_TYPE_CASTER(akantu::ArrayProxy<akantu::Real>, _("ArrayProxy"));
//     bool load(handle, bool) {
//       // /* Extract PyObject from handle */
//       // PyObject *source = src.ptr();
//       // /* Try converting into a Python integer value */
//       // PyObject *tmp = PyNumber_Long(source);
//       // if (!tmp)
//       //   return false;
//       // /* Now try to convert into a C++ int */
//       // value.long_value = PyLong_AsLong(tmp);
//       // Py_DECREF(tmp);
//       // /* Ensure return code was OK (to avoid out-of-range errors etc) */
//       // return !(value.long_value == -1 && !PyErr_Occurred());
//       return false;
//     }

//     static handle cast(akantu::ArrayProxy<akantu::Real> & src,
//                        return_value_policy /* policy */, handle parent) {
//       constexpr ssize_t elem_size = sizeof(akantu::Real);
//       ssize_t nb_comp = src.getNbComponent();
//       ssize_t size = src.size();
//       std::cout << "<memory at " << reinterpret_cast<void *>(src.storage())
//                 << ">\n";

//       auto a = array({size, nb_comp}, {elem_size * nb_comp, elem_size},
//                      src.storage(), handle());
//       return a.release();
//     }
//   };
// } // namespace detail
// } //  namespace pybind11

#endif /* __AKANTU_PYBIND11_AKANTU_HH__ */
