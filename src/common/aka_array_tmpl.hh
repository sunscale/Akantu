/**
 * @file   aka_array_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Jul 15 2010
 * @date last modification: Fri Jan 22 2016
 *
 * @brief  Inline functions of the classes Array<T> and ArrayBase
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
/* Inline Functions Array<T>                                                 */
/* -------------------------------------------------------------------------- */
#include "aka_array.hh"
/* -------------------------------------------------------------------------- */
#include <memory>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_AKA_ARRAY_TMPL_HH__
#define __AKANTU_AKA_ARRAY_TMPL_HH__

namespace akantu {

namespace debug {
  struct ArrayException : public Exception {};
} // namespace debug

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline T & Array<T, is_scal>::operator()(UInt i, UInt j) {
  AKANTU_DEBUG_ASSERT(size_ > 0, "The array \"" << id << "\" is empty");
  AKANTU_DEBUG_ASSERT((i < size_) && (j < nb_component),
                      "The value at position ["
                          << i << "," << j << "] is out of range in array \""
                          << id << "\"");
  return values[i * nb_component + j];
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline const T & Array<T, is_scal>::operator()(UInt i, UInt j) const {
  AKANTU_DEBUG_ASSERT(size_ > 0, "The array \"" << id << "\" is empty");
  AKANTU_DEBUG_ASSERT((i < size_) && (j < nb_component),
                      "The value at position ["
                          << i << "," << j << "] is out of range in array \""
                          << id << "\"");
  return values[i * nb_component + j];
}

template <class T, bool is_scal>
inline T & Array<T, is_scal>::operator[](UInt i) {
  AKANTU_DEBUG_ASSERT(size_ > 0, "The array \"" << id << "\" is empty");
  AKANTU_DEBUG_ASSERT((i < size_ * nb_component),
                      "The value at position ["
                          << i << "] is out of range in array \"" << id
                          << "\"");
  return values[i];
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline const T & Array<T, is_scal>::operator[](UInt i) const {
  AKANTU_DEBUG_ASSERT(size_ > 0, "The array \"" << id << "\" is empty");
  AKANTU_DEBUG_ASSERT((i < size_ * nb_component),
                      "The value at position ["
                          << i << "] is out of range in array \"" << id
                          << "\"");
  return values[i];
}

/* -------------------------------------------------------------------------- */
/**
 * append a tuple to the array with the value value for all components
 * @param value the new last tuple or the array will contain nb_component copies
 * of value
 */
template <class T, bool is_scal>
inline void Array<T, is_scal>::push_back(const T & value) {
  resizeUnitialized(size_ + 1, true, value);
}

/* -------------------------------------------------------------------------- */
/**
 * append a tuple to the array
 * @param new_elem a C-array containing the values to be copied to the end of
 * the array */
// template <class T, bool is_scal>
// inline void Array<T, is_scal>::push_back(const T new_elem[]) {
//   UInt pos = size_;

//   resizeUnitialized(size_ + 1, false);

//   T * tmp = values + nb_component * pos;
//   std::uninitialized_copy(new_elem, new_elem + nb_component, tmp);
// }

/* -------------------------------------------------------------------------- */
#ifndef SWIG
/**
 * append a matrix or a vector to the array
 * @param new_elem a reference to a Matrix<T> or Vector<T> */
template <class T, bool is_scal>
template <template <typename> class C, typename>
inline void Array<T, is_scal>::push_back(const C<T> & new_elem) {
  AKANTU_DEBUG_ASSERT(
      nb_component == new_elem.size(),
      "The vector("
          << new_elem.size()
          << ") as not a size compatible with the Array (nb_component="
          << nb_component << ").");
  UInt pos = size_;
  resizeUnitialized(size_ + 1, false);

  T * tmp = values + nb_component * pos;
  std::uninitialized_copy(new_elem.storage(), new_elem.storage() + nb_component,
                          tmp);
}

/* -------------------------------------------------------------------------- */
/**
 * append a tuple to the array
 * @param it an iterator to the tuple to be copied to the end of the array */
template <class T, bool is_scal>
template <class Ret>
inline void
Array<T, is_scal>::push_back(const Array<T, is_scal>::iterator<Ret> & it) {
  UInt pos = size_;

  resizeUnitialized(size_ + 1, false);

  T * tmp = values + nb_component * pos;
  T * new_elem = it.data();
  std::uninitialized_copy(new_elem, new_elem + nb_component, tmp);
}
#endif
/* -------------------------------------------------------------------------- */
/**
 * erase an element. If the erased element is not the last of the array, the
 * last element is moved into the hole in order to maintain contiguity. This
 * may invalidate existing iterators (For instance an iterator obtained by
 * Array::end() is no longer correct) and will change the order of the
 * elements.
 * @param i index of element to erase
 */
template <class T, bool is_scal> inline void Array<T, is_scal>::erase(UInt i) {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT((size_ > 0), "The array is empty");
  AKANTU_DEBUG_ASSERT((i < size_),
                      "The element at position [" << i << "] is out of range ("
                                                  << i << ">=" << size_ << ")");

  if (i != (size_ - 1)) {
    for (UInt j = 0; j < nb_component; ++j) {
      values[i * nb_component + j] = values[(size_ - 1) * nb_component + j];
    }
  }

  resize(size_ - 1);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Subtract another array entry by entry from this array in place. Both arrays
 * must
 * have the same size and nb_component. If the arrays have different shapes,
 * code compiled in debug mode will throw an expeption and optimised code
 * will behave in an unpredicted manner
 * @param other array to subtract from this
 * @return reference to modified this
 */
template <class T, bool is_scal>
Array<T, is_scal> & Array<T, is_scal>::
operator-=(const Array<T, is_scal> & vect) {
  AKANTU_DEBUG_ASSERT((size_ == vect.size_) &&
                          (nb_component == vect.nb_component),
                      "The too array don't have the same sizes");

  T * a = values;
  T * b = vect.storage();
  for (UInt i = 0; i < size_ * nb_component; ++i) {
    *a -= *b;
    ++a;
    ++b;
  }

  return *this;
}

/* -------------------------------------------------------------------------- */
/**
 * Add another array entry by entry to this array in place. Both arrays must
 * have the same size and nb_component. If the arrays have different shapes,
 * code compiled in debug mode will throw an expeption and optimised code
 * will behave in an unpredicted manner
 * @param other array to add to this
 * @return reference to modified this
 */
template <class T, bool is_scal>
Array<T, is_scal> & Array<T, is_scal>::
operator+=(const Array<T, is_scal> & vect) {
  AKANTU_DEBUG_ASSERT((size_ == vect.size) &&
                          (nb_component == vect.nb_component),
                      "The too array don't have the same sizes");

  T * a = values;
  T * b = vect.storage();
  for (UInt i = 0; i < size_ * nb_component; ++i) {
    *a++ += *b++;
  }

  return *this;
}

/* -------------------------------------------------------------------------- */
/**
 * Multiply all entries of this array by a scalar in place
 * @param alpha scalar multiplicant
 * @return reference to modified this
 */
template <class T, bool is_scal>
Array<T, is_scal> & Array<T, is_scal>::operator*=(const T & alpha) {
  T * a = values;
  for (UInt i = 0; i < size_ * nb_component; ++i) {
    *a++ *= alpha;
  }

  return *this;
}

/* -------------------------------------------------------------------------- */
/**
 * Compare this array element by element to another.
 * @param other array to compare to
 * @return true it all element are equal and arrays have the same shape, else
 * false
 */
template <class T, bool is_scal>
bool Array<T, is_scal>::operator==(const Array<T, is_scal> & array) const {
  bool equal = nb_component == array.nb_component && size_ == array.size_ &&
               id == array.id;
  if (!equal)
    return false;

  if (values == array.storage())
    return true;
  else
    return std::equal(values, values + size_ * nb_component, array.storage());
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
bool Array<T, is_scal>::operator!=(const Array<T, is_scal> & array) const {
  return !operator==(array);
}

/* -------------------------------------------------------------------------- */
#ifndef SWIG
/**
 * set all tuples of the array to a given vector or matrix
 * @param vm Matrix or Vector to fill the array with
 */
template <class T, bool is_scal>
template <template <typename> class C, typename>
inline void Array<T, is_scal>::set(const C<T> & vm) {
  AKANTU_DEBUG_ASSERT(
      nb_component == vm.size(),
      "The size of the object does not match the number of components");
  for (T * it = values; it < values + nb_component * size_;
       it += nb_component) {
    std::copy(vm.storage(), vm.storage() + nb_component, it);
  }
}
#endif
/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
void Array<T, is_scal>::append(const Array<T> & other) {
  AKANTU_DEBUG_ASSERT(
      nb_component == other.nb_component,
      "Cannot append an array with a different number of component");
  UInt old_size = this->size_;
  this->resizeUnitialized(this->size_ + other.size(), false);

  T * tmp = values + nb_component * old_size;
  std::uninitialized_copy(other.storage(),
                          other.storage() + other.size() * nb_component, tmp);
}

/* -------------------------------------------------------------------------- */
/* Functions Array<T, is_scal>                                                */
/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
Array<T, is_scal>::Array(UInt size, UInt nb_component, const ID & id)
    : ArrayBase(id), values(nullptr) {
  AKANTU_DEBUG_IN();
  allocate(size, nb_component);

  if (!is_scal) {
    T val = T();
    std::uninitialized_fill(values, values + size * nb_component, val);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
Array<T, is_scal>::Array(UInt size, UInt nb_component, const T def_values[],
                         const ID & id)
    : ArrayBase(id), values(NULL) {
  AKANTU_DEBUG_IN();
  allocate(size, nb_component);

  T * tmp = values;

  for (UInt i = 0; i < size; ++i) {
    tmp = values + nb_component * i;
    std::uninitialized_copy(def_values, def_values + nb_component, tmp);
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
Array<T, is_scal>::Array(UInt size, UInt nb_component, const T & value,
                         const ID & id)
    : ArrayBase(id), values(nullptr) {
  AKANTU_DEBUG_IN();
  allocate(size, nb_component);

  std::uninitialized_fill_n(values, size * nb_component, value);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
Array<T, is_scal>::Array(const Array<T, is_scal> & vect, bool deep,
                         const ID & id)
    : ArrayBase(vect) {
  AKANTU_DEBUG_IN();
  this->id = (id == "") ? vect.id : id;

  if (deep) {
    allocate(vect.size_, vect.nb_component);
    T * tmp = values;
    std::uninitialized_copy(vect.storage(),
                            vect.storage() + size_ * nb_component, tmp);
  } else {
    this->values = vect.storage();
    this->size_ = vect.size_;
    this->nb_component = vect.nb_component;
    this->allocated_size = vect.allocated_size;
    this->size_of_type = vect.size_of_type;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
#ifndef SWIG
template <class T, bool is_scal>
Array<T, is_scal>::Array(const std::vector<T> & vect) {
  AKANTU_DEBUG_IN();
  this->id = "";

  allocate(vect.size(), 1);
  T * tmp = values;
  std::uninitialized_copy(&(vect[0]), &(vect[size_ - 1]), tmp);

  AKANTU_DEBUG_OUT();
}
#endif

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal> Array<T, is_scal>::~Array() {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG(dblAccessory,
               "Freeing " << printMemorySize<T>(allocated_size * nb_component)
                          << " (" << id << ")");

  if (values) {
    if (!is_scal)
      for (UInt i = 0; i < size_ * nb_component; ++i) {
        T * obj = values + i;
        obj->~T();
      }
    free(values);
  }
  size_ = allocated_size = 0;
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * perform the allocation for the constructors
 * @param size is the size of the array
 * @param nb_component is the number of component of the array
 */
template <class T, bool is_scal>
void Array<T, is_scal>::allocate(UInt size, UInt nb_component) {
  AKANTU_DEBUG_IN();
  if (size == 0) {
    values = nullptr;
  } else {
    values = static_cast<T *>(malloc(nb_component * size * sizeof(T)));
    AKANTU_DEBUG_ASSERT(values != nullptr,
                        "Cannot allocate "
                            << printMemorySize<T>(size * nb_component) << " ("
                            << id << ")");
  }

  if (values == nullptr) {
    this->size_ = this->allocated_size = 0;
  } else {
    AKANTU_DEBUG(dblAccessory,
                 "Allocated " << printMemorySize<T>(size * nb_component) << " ("
                              << id << ")");
    this->size_ = this->allocated_size = size;
  }

  this->size_of_type = sizeof(T);
  this->nb_component = nb_component;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
void Array<T, is_scal>::reserve(UInt new_size) {
  UInt tmp_size = this->size_;
  resizeUnitialized(new_size, false);
  this->size_ = tmp_size;
}

/* -------------------------------------------------------------------------- */
/**
 * change the size of the array and allocate or free memory if needed. If the
 * size increases, the new tuples are filled with zeros
 * @param new_size new number of tuples contained in the array */
template <class T, bool is_scal> void Array<T, is_scal>::resize(UInt new_size) {
  resizeUnitialized(new_size, !is_scal);
}

/* -------------------------------------------------------------------------- */
/**
 * change the size of the array and allocate or free memory if needed. If the
 * size increases, the new tuples are filled with zeros
 * @param new_size new number of tuples contained in the array */
template <class T, bool is_scal>
void Array<T, is_scal>::resize(UInt new_size, const T & val) {
  this->resizeUnitialized(new_size, true, val);
}

/* -------------------------------------------------------------------------- */
/**
 * change the size of the array and allocate or free memory if needed.
 * @param new_size new number of tuples contained in the array */
template <class T, bool is_scal>
void Array<T, is_scal>::resizeUnitialized(UInt new_size, bool fill,
                                          const T & val) {
  //  AKANTU_DEBUG_IN();
  // free some memory
  if (new_size <= allocated_size) {
    if (!is_scal) {
      T * old_values = values;
      if (new_size < size_) {
        for (UInt i = new_size * nb_component; i < size_ * nb_component; ++i) {
          T * obj = old_values + i;
          obj->~T();
        }
      }
    }

    if (allocated_size - new_size > AKANTU_MIN_ALLOCATION) {
      AKANTU_DEBUG(dblAccessory,
                   "Freeing " << printMemorySize<T>((allocated_size - size_) *
                                                    nb_component)
                              << " (" << id << ")");

      // Normally there are no allocation problem when reducing an array
      if (new_size == 0) {
        free(values);
        values = nullptr;
      } else {
        auto * tmp_ptr = static_cast<T *>(
            realloc(values, new_size * nb_component * sizeof(T)));

        if (tmp_ptr == nullptr) {
          AKANTU_EXCEPTION("Cannot free data ("
                           << id << ")"
                           << " [current allocated size : " << allocated_size
                           << " | "
                           << "requested size : " << new_size << "]");
        }
        values = tmp_ptr;
      }
      allocated_size = new_size;
    }
  } else {
    // allocate more memory
    UInt size_to_alloc = (new_size - allocated_size < AKANTU_MIN_ALLOCATION)
                             ? allocated_size + AKANTU_MIN_ALLOCATION
                             : new_size;

    auto * tmp_ptr = static_cast<T *>(
        realloc(values, size_to_alloc * nb_component * sizeof(T)));
    AKANTU_DEBUG_ASSERT(
        tmp_ptr != nullptr,
        "Cannot allocate " << printMemorySize<T>(size_to_alloc * nb_component));
    if (tmp_ptr == nullptr) {
      AKANTU_DEBUG_ERROR("Cannot allocate more data ("
                         << id << ")"
                         << " [current allocated size : " << allocated_size
                         << " | "
                         << "requested size : " << new_size << "]");
    }

    AKANTU_DEBUG(dblAccessory,
                 "Allocating " << printMemorySize<T>(
                     (size_to_alloc - allocated_size) * nb_component));

    allocated_size = size_to_alloc;
    values = tmp_ptr;
  }

  if (fill && this->size_ < new_size) {
    std::uninitialized_fill(values + size_ * nb_component,
                            values + new_size * nb_component, val);
  }

  size_ = new_size;
  //  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * change the number of components by interlacing data
 * @param multiplicator number of interlaced components add
 * @param block_size blocks of data in the array
 * Examaple for block_size = 2, multiplicator = 2
 * array = oo oo oo -> new array = oo nn nn oo nn nn oo nn nn */
template <class T, bool is_scal>
void Array<T, is_scal>::extendComponentsInterlaced(UInt multiplicator,
                                                   UInt block_size) {
  AKANTU_DEBUG_IN();

  if (multiplicator == 1)
    return;

  AKANTU_DEBUG_ASSERT(multiplicator > 1, "invalid multiplicator");
  AKANTU_DEBUG_ASSERT(nb_component % block_size == 0,
                      "stride must divide actual number of components");

  values = static_cast<T *>(
      realloc(values, nb_component * multiplicator * size_ * sizeof(T)));

  UInt new_component = nb_component / block_size * multiplicator;

  for (UInt i = 0, k = size_ - 1; i < size_; ++i, --k) {
    for (UInt j = 0; j < new_component; ++j) {
      UInt m = new_component - j - 1;
      UInt n = m / multiplicator;
      for (UInt l = 0, p = block_size - 1; l < block_size; ++l, --p) {
        values[k * nb_component * multiplicator + m * block_size + p] =
            values[k * nb_component + n * block_size + p];
      }
    }
  }

  nb_component = nb_component * multiplicator;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * search elem in the array, return  the position of the first occurrence or
 * -1 if not found
 *  @param elem the element to look for
 *  @return index of the first occurrence of elem or -1 if elem is not present
 */
template <class T, bool is_scal>
UInt Array<T, is_scal>::find(const T & elem) const {
  AKANTU_DEBUG_IN();

  auto begin = this->begin();
  auto end = this->end();
  auto it = std::find(begin, end, elem);

  AKANTU_DEBUG_OUT();
  return (it != end) ? it - begin : UInt(-1);
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal> UInt Array<T, is_scal>::find(T elem[]) const {
  AKANTU_DEBUG_IN();
  T * it = values;
  UInt i = 0;
  for (; i < size_; ++i) {
    if (*it == elem[0]) {
      T * cit = it;
      UInt c = 0;
      for (; (c < nb_component) && (*cit == elem[c]); ++c, ++cit)
        ;
      if (c == nb_component) {
        AKANTU_DEBUG_OUT();
        return i;
      }
    }
    it += nb_component;
  }
  return UInt(-1);
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
template <template <typename> class C, typename>
inline UInt Array<T, is_scal>::find(const C<T> & elem) {
  AKANTU_DEBUG_ASSERT(elem.size() == nb_component,
                      "Cannot find an element with a wrong size ("
                          << elem.size() << ") != " << nb_component);
  return this->find(elem.storage());
}

/* -------------------------------------------------------------------------- */
/**
 * copy the content of another array. This overwrites the current content.
 * @param other Array to copy into this array. It has to have the same
 * nb_component as this. If compiled in debug mode, an incorrect other will
 * result in an exception being thrown. Optimised code may result in
 * unpredicted behaviour.
 */
template <class T, bool is_scal>
void Array<T, is_scal>::copy(const Array<T, is_scal> & vect,
                             bool no_sanity_check) {
  AKANTU_DEBUG_IN();

  if (!no_sanity_check)
    if (vect.nb_component != nb_component)
      AKANTU_DEBUG_ERROR(
          "The two arrays do not have the same number of components");

  resize((vect.size_ * vect.nb_component) / nb_component);

  T * tmp = values;
  std::uninitialized_copy(vect.storage(), vect.storage() + size_ * nb_component,
                          tmp);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <bool is_scal> class ArrayPrintHelper {
public:
  template <typename T>
  static void print_content(const Array<T> & vect, std::ostream & stream,
                            int indent) {
    if (AKANTU_DEBUG_TEST(dblDump) || AKANTU_DEBUG_LEVEL_IS_TEST()) {
      std::string space;
      for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
        ;

      stream << space << " + values         : {";
      for (UInt i = 0; i < vect.size(); ++i) {
        stream << "{";
        for (UInt j = 0; j < vect.getNbComponent(); ++j) {
          stream << vect(i, j);
          if (j != vect.getNbComponent() - 1)
            stream << ", ";
        }
        stream << "}";
        if (i != vect.size() - 1)
          stream << ", ";
      }
      stream << "}" << std::endl;
    }
  }
};

template <> class ArrayPrintHelper<false> {
public:
  template <typename T>
  static void print_content(__attribute__((unused)) const Array<T> & vect,
                            __attribute__((unused)) std::ostream & stream,
                            __attribute__((unused)) int indent) {}
};

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
void Array<T, is_scal>::printself(std::ostream & stream, int indent) const {
  std::string space;
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
    ;

  std::streamsize prec = stream.precision();
  std::ios_base::fmtflags ff = stream.flags();

  stream.setf(std::ios_base::showbase);
  stream.precision(2);

  stream << space << "Array<" << debug::demangle(typeid(T).name()) << "> ["
         << std::endl;
  stream << space << " + id             : " << this->id << std::endl;
  stream << space << " + size           : " << this->size_ << std::endl;
  stream << space << " + nb_component   : " << this->nb_component << std::endl;
  stream << space << " + allocated size : " << this->allocated_size
         << std::endl;
  stream << space << " + memory size    : "
         << printMemorySize<T>(allocated_size * nb_component) << std::endl;
  if (!AKANTU_DEBUG_LEVEL_IS_TEST())
    stream << space << " + address        : " << std::hex << this->values
           << std::dec << std::endl;

  stream.precision(prec);
  stream.flags(ff);

  ArrayPrintHelper<is_scal>::print_content(*this, stream, indent);

  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
/* Inline Functions ArrayBase                                                */
/* -------------------------------------------------------------------------- */

inline UInt ArrayBase::getMemorySize() const {
  return allocated_size * nb_component * size_of_type;
}

inline void ArrayBase::empty() { size_ = 0; }

#ifndef SWIG
/* -------------------------------------------------------------------------- */
/* Iterators                                                                  */
/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
template <class R, class daughter, class IR, bool is_tensor>
class Array<T, is_scal>::iterator_internal {
public:
  using value_type = R;
  using pointer = R *;
  using reference = R &;
  using proxy = typename R::proxy;
  using const_proxy = const typename R::proxy;
  using const_reference = const R &;
  using internal_value_type = IR;
  using internal_pointer = IR *;
  using difference_type = std::ptrdiff_t;
  using iterator_category = std::random_access_iterator_tag;

public:
  iterator_internal() = default;

  iterator_internal(pointer_type data, UInt _offset)
      : _offset(_offset), initial(data), ret(nullptr), ret_ptr(data) {
    AKANTU_DEBUG_ERROR(
        "The constructor should never be called it is just an ugly trick...");
  }

  iterator_internal(std::unique_ptr<internal_value_type> && wrapped)
      : _offset(wrapped->size()), initial(wrapped->storage()),
        ret(std::move(wrapped)), ret_ptr(ret->storage()) {}

  iterator_internal(const iterator_internal & it) {
    if (this != &it) {
      this->_offset = it._offset;
      this->initial = it.initial;
      this->ret_ptr = it.ret_ptr;
      this->ret = std::make_unique<internal_value_type>(*it.ret, false);
    }
  }

  iterator_internal(iterator_internal && it) = default;

  virtual ~iterator_internal() = default;

  inline iterator_internal & operator=(const iterator_internal & it) {
    if (this != &it) {
      this->_offset = it._offset;
      this->initial = it.initial;
      this->ret_ptr = it.ret_ptr;
      if (this->ret)
        this->ret->shallowCopy(*it.ret);
      else
        this->ret = std::make_unique<internal_value_type>(*it.ret, false);
    }
    return *this;
  }

  UInt getCurrentIndex() {
    return (this->ret_ptr - this->initial) / this->_offset;
  };

  inline reference operator*() {
    ret->values = ret_ptr;
    return *ret;
  };
  inline const_reference operator*() const {
    ret->values = ret_ptr;
    return *ret;
  };
  inline pointer operator->() {
    ret->values = ret_ptr;
    return ret.get();
  };
  inline daughter & operator++() {
    ret_ptr += _offset;
    return static_cast<daughter &>(*this);
  };
  inline daughter & operator--() {
    ret_ptr -= _offset;
    return static_cast<daughter &>(*this);
  };

  inline daughter & operator+=(const UInt n) {
    ret_ptr += _offset * n;
    return static_cast<daughter &>(*this);
  }
  inline daughter & operator-=(const UInt n) {
    ret_ptr -= _offset * n;
    return static_cast<daughter &>(*this);
  }

  inline proxy operator[](const UInt n) {
    ret->values = ret_ptr + n * _offset;
    return proxy(*ret);
  }
  inline const_proxy operator[](const UInt n) const {
    ret->values = ret_ptr + n * _offset;
    return const_proxy(*ret);
  }

  inline bool operator==(const iterator_internal & other) const {
    return this->ret_ptr == other.ret_ptr;
  }
  inline bool operator!=(const iterator_internal & other) const {
    return this->ret_ptr != other.ret_ptr;
  }
  inline bool operator<(const iterator_internal & other) const {
    return this->ret_ptr < other.ret_ptr;
  }
  inline bool operator<=(const iterator_internal & other) const {
    return this->ret_ptr <= other.ret_ptr;
  }
  inline bool operator>(const iterator_internal & other) const {
    return this->ret_ptr > other.ret_ptr;
  }
  inline bool operator>=(const iterator_internal & other) const {
    return this->ret_ptr >= other.ret_ptr;
  }

  inline daughter operator+(difference_type n) {
    daughter tmp(static_cast<daughter &>(*this));
    tmp += n;
    return tmp;
  }
  inline daughter operator-(difference_type n) {
    daughter tmp(static_cast<daughter &>(*this));
    tmp -= n;
    return tmp;
  }

  inline difference_type operator-(const iterator_internal & b) {
    return (this->ret_ptr - b.ret_ptr) / _offset;
  }

  inline pointer_type data() const { return ret_ptr; }
  inline difference_type offset() const { return _offset; }

protected:
  UInt _offset{0};
  pointer_type initial{nullptr};
  std::unique_ptr<internal_value_type> ret{nullptr};
  pointer_type ret_ptr{nullptr};
};

/* -------------------------------------------------------------------------- */
/**
 * Specialization for scalar types
 */
template <class T, bool is_scal>
template <class R, class daughter, class IR>
class Array<T, is_scal>::iterator_internal<R, daughter, IR, false> {
public:
  using value_type = R;
  using pointer = R *;
  using reference = R &;
  using const_reference = const R &;
  using internal_value_type = IR;
  using internal_pointer = IR *;
  using difference_type = std::ptrdiff_t;
  using iterator_category = std::random_access_iterator_tag;

public:
  iterator_internal(pointer data = nullptr) : ret(data), initial(data){};
  iterator_internal(const iterator_internal & it) = default;
  iterator_internal(iterator_internal && it) = default;

  virtual ~iterator_internal() = default;

  inline iterator_internal & operator=(const iterator_internal & it) = default;

  UInt getCurrentIndex() { return (this->ret - this->initial); };

  inline reference operator*() { return *ret; };
  inline const_reference operator*() const { return *ret; };
  inline pointer operator->() { return ret; };
  inline daughter & operator++() {
    ++ret;
    return static_cast<daughter &>(*this);
  };
  inline daughter & operator--() {
    --ret;
    return static_cast<daughter &>(*this);
  };

  inline daughter & operator+=(const UInt n) {
    ret += n;
    return static_cast<daughter &>(*this);
  }
  inline daughter & operator-=(const UInt n) {
    ret -= n;
    return static_cast<daughter &>(*this);
  }

  inline reference operator[](const UInt n) { return ret[n]; }

  inline bool operator==(const iterator_internal & other) const {
    return ret == other.ret;
  }
  inline bool operator!=(const iterator_internal & other) const {
    return ret != other.ret;
  }
  inline bool operator<(const iterator_internal & other) const {
    return ret < other.ret;
  }
  inline bool operator<=(const iterator_internal & other) const {
    return ret <= other.ret;
  }
  inline bool operator>(const iterator_internal & other) const {
    return ret > other.ret;
  }
  inline bool operator>=(const iterator_internal & other) const {
    return ret >= other.ret;
  }

  inline daughter operator-(difference_type n) { return daughter(ret - n); }
  inline daughter operator+(difference_type n) { return daughter(ret + n); }

  inline difference_type operator-(const iterator_internal & b) {
    return ret - b.ret;
  }

  inline pointer data() const { return ret; }

protected:
  pointer ret{nullptr};
  pointer initial{nullptr};
};

/* -------------------------------------------------------------------------- */
/* Begin/End functions implementation                                         */
/* -------------------------------------------------------------------------- */
namespace detail {
  template <class Tuple, size_t... Is>
  constexpr auto take_front_impl(Tuple && t, std::index_sequence<Is...>) {
    return std::make_tuple(std::get<Is>(std::forward<Tuple>(t))...);
  }

  template <size_t N, class Tuple> constexpr auto take_front(Tuple && t) {
    return take_front_impl(std::forward<Tuple>(t),
                           std::make_index_sequence<N>{});
  }

  template <typename... V>
  constexpr auto product_all(V &&... v) ->
      typename std::common_type<V...>::type {
    typename std::common_type<V...>::type result = 1;
    (void)std::initializer_list<int>{(result *= v, 0)...};
    return result;
  }

  template <typename... T> std::string to_string_all(T &&... t) {
    if (sizeof...(T) == 0)
      return "";

    std::stringstream ss;
    bool noComma = true;
    ss << "(";
    (void)std::initializer_list<bool>{
        (ss << (noComma ? "" : ", ") << t, noComma = false)...};
    ss << ")";
    return ss.str();
  }

  template <std::size_t N> struct InstantiationHelper {
    template <typename type, typename T, typename... Ns>
    static auto instantiate(T && data, Ns... ns) {
      return std::make_unique<type>(data, ns...);
    }
  };

  template <> struct InstantiationHelper<0> {
    template <typename type, typename T> static auto instantiate(T && data) {
      return data;
    }
  };

  template <typename Arr, typename T, typename... Ns>
  decltype(auto) get_iterator(Arr && array, T * data, Ns &&... ns) {
    using type = IteratorHelper_t<sizeof...(Ns) - 1, T>;
    using array_type = std::decay_t<Arr>;
    using iterator =
        std::conditional_t<std::is_const<std::remove_reference_t<Arr>>::value,
                           typename array_type::template const_iterator<type>,
                           typename array_type::template iterator<type>>;
    static_assert(sizeof...(Ns), "You should provide a least one size");

    if (array.getNbComponent() * array.size() !=
        product_all(std::forward<Ns>(ns)...)) {
      AKANTU_CUSTOM_EXCEPTION_INFO(
          debug::ArrayException(),
          "The iterator on " << debug::demangle(typeid(Arr).name())
              << to_string_all(array.size(), array.getNbComponent())
              << "is not compatible with the type "
              << debug::demangle(typeid(type).name()) << to_string_all(ns...));
    }

    auto && wrapped = aka::apply(
        [&](auto... n) {
          return InstantiationHelper<sizeof...(n)>::template instantiate<type>(
              data, n...);
        },
        take_front<sizeof...(Ns) - 1>(std::make_tuple(ns...)));

    return iterator(std::move(wrapped));
  }
} // namespace detail

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
template <typename... Ns>
inline decltype(auto) Array<T, is_scal>::begin(Ns &&... ns) {
  return detail::get_iterator(*this, values, std::forward<Ns>(ns)..., size_);
}

template <class T, bool is_scal>
template <typename... Ns>
inline decltype(auto) Array<T, is_scal>::end(Ns &&... ns) {
  return detail::get_iterator(*this, values + nb_component * size_,
                              std::forward<Ns>(ns)..., size_);
}

template <class T, bool is_scal>
template <typename... Ns>
inline decltype(auto) Array<T, is_scal>::begin(Ns &&... ns) const {
  return detail::get_iterator(*this, values, std::forward<Ns>(ns)..., size_);
}

template <class T, bool is_scal>
template <typename... Ns>
inline decltype(auto) Array<T, is_scal>::end(Ns &&... ns) const {
  return detail::get_iterator(*this, values + nb_component * size_,
                              std::forward<Ns>(ns)..., size_);
}

template <class T, bool is_scal>
template <typename... Ns>
inline decltype(auto) Array<T, is_scal>::begin_reinterpret(Ns &&... ns) {
  return detail::get_iterator(*this, values, std::forward<Ns>(ns)...);
}

template <class T, bool is_scal>
template <typename... Ns>
inline decltype(auto) Array<T, is_scal>::end_reinterpret(Ns &&... ns) {
  return detail::get_iterator(
      *this, values + detail::product_all(std::forward<Ns>(ns)...),
      std::forward<Ns>(ns)...);
}

template <class T, bool is_scal>
template <typename... Ns>
inline decltype(auto) Array<T, is_scal>::begin_reinterpret(Ns &&... ns) const {
  return detail::get_iterator(*this, values, std::forward<Ns>(ns)...);
}

template <class T, bool is_scal>
template <typename... Ns>
inline decltype(auto) Array<T, is_scal>::end_reinterpret(Ns &&... ns) const {
  return detail::get_iterator(
      *this, values + detail::product_all(std::forward<Ns>(ns)...),
      std::forward<Ns>(ns)...);
}

/* -------------------------------------------------------------------------- */
/* Views                                                                      */
/* -------------------------------------------------------------------------- */
namespace detail {
  template <typename Array, typename... Ns> class ArrayView {
    using tuple = std::tuple<Ns...>;

  public:
    ArrayView(Array && array, Ns... ns)
        : array(std::forward<Array>(array)), sizes(std::move(ns)...){};

    decltype(auto) begin() {
      return aka::apply(
          [&](auto &&... ns) { return array.begin_reinterpret(ns...); }, sizes);
    }

    decltype(auto) end() {
      return aka::apply(
          [&](auto &&... ns) { return array.end_reinterpret(ns...); }, sizes);
    }

    decltype(auto) size() {
      return std::get<std::tuple_size<tuple>::value - 1>(sizes);
    }

    decltype(auto) dims() {
      return std::tuple_size<tuple>::value - 1;
    }
  private:
    Array array;
    tuple sizes;
  };
} // namespace detail

/* -------------------------------------------------------------------------- */
template <typename Array, typename... Ns>
decltype(auto) make_view(Array && array, Ns... ns) {
  static_assert(sizeof...(Ns), "You should provide a least one dimension");
  auto size = std::forward<decltype(array)>(array).size() *
              std::forward<decltype(array)>(array).getNbComponent() /
              detail::product_all(ns...);

  return detail::ArrayView<Array, std::common_type_t<size_t, Ns>...,
                           std::common_type_t<size_t, decltype(size)>>(
      std::forward<Array>(array), std::move(ns)..., size);
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
template <typename R>
class Array<T, is_scal>::const_iterator
    : public iterator_internal<const R, Array<T, is_scal>::const_iterator<R>,
                               R> {
public:
  using parent = iterator_internal<const R, const_iterator, R>;
  using value_type = typename parent::value_type;
  using pointer = typename parent::pointer;
  using reference = typename parent::reference;
  using difference_type = typename parent::difference_type;
  using iterator_category = typename parent::iterator_category;

public:
  const_iterator() : parent(){};
  // const_iterator(pointer_type data, UInt offset) : parent(data, offset) {}
  // const_iterator(pointer warped) : parent(warped) {}
  // const_iterator(const parent & it) : parent(it) {}

  const_iterator(const const_iterator & it) = default;
  const_iterator(const_iterator && it) = default;

  template <typename P, typename = std::enable_if_t<not is_tensor<P>::value>>
  const_iterator(P * data) : parent(data) {}

  template <typename UP_P,
            typename =
                std::enable_if_t<is_tensor<typename UP_P::element_type>::value>>
  const_iterator(UP_P && tensor) : parent(std::forward<UP_P>(tensor)) {}

  const_iterator & operator=(const const_iterator & it) = default;
};

template <class T, class R, bool is_tensor_ = is_tensor<R>::value>
struct ConstConverterIteratorHelper {
  using const_iterator = typename Array<T>::template const_iterator<R>;
  using iterator = typename Array<T>::template iterator<R>;

  static inline const_iterator convert(const iterator & it) {
    return const_iterator(std::unique_ptr<R>(new R(*it, false)));
  }
};

template <class T, class R> struct ConstConverterIteratorHelper<T, R, false> {
  using const_iterator = typename Array<T>::template const_iterator<R>;
  using iterator = typename Array<T>::template iterator<R>;
  static inline const_iterator convert(const iterator & it) {
    return const_iterator(it.data());
  }
};

template <class T, bool is_scal>
template <typename R>
class Array<T, is_scal>::iterator
    : public iterator_internal<R, Array<T, is_scal>::iterator<R>> {
public:
  using parent = iterator_internal<R, iterator>;
  using value_type = typename parent::value_type;
  using pointer = typename parent::pointer;
  using reference = typename parent::reference;
  using difference_type = typename parent::difference_type;
  using iterator_category = typename parent::iterator_category;

public:
  iterator() : parent(){};
  iterator(const iterator & it) = default;
  iterator(iterator && it) = default;

  template <typename P, typename = std::enable_if_t<not is_tensor<P>::value>>
  iterator(P * data) : parent(data) {}

  template <typename UP_P,
            typename =
                std::enable_if_t<is_tensor<typename UP_P::element_type>::value>>
  iterator(UP_P && tensor) : parent(std::forward<UP_P>(tensor)) {}

  iterator & operator=(const iterator & it) = default;

  operator const_iterator<R>() {
    return ConstConverterIteratorHelper<T, R>::convert(*this);
  }
};

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
template <typename R>
inline typename Array<T, is_scal>::template iterator<R>
Array<T, is_scal>::erase(const iterator<R> & it) {
  T * curr = it.data();
  UInt pos = (curr - values) / nb_component;
  erase(pos);
  iterator<R> rit = it;
  return --rit;
}
#endif

} // namespace akantu

#endif /* __AKANTU_AKA_ARRAY_TMPL_HH__ */
