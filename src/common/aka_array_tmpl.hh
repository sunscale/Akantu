/**
 * @file   aka_array_tmpl.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Jul 15 2010
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Inline functions of the classes Array<T> and ArrayBase
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
/* Inline Functions Array<T>                                                  */
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
/* -------------------------------------------------------------------------- */
template <typename T, ArrayAllocationType allocation_trait>
ArrayDataLayer<T, allocation_trait>::ArrayDataLayer(UInt size,
                                                    UInt nb_component,
                                                    const ID & id)
    : ArrayBase(id) {
  allocate(size, nb_component);
}

/* -------------------------------------------------------------------------- */
template <typename T, ArrayAllocationType allocation_trait>
ArrayDataLayer<T, allocation_trait>::ArrayDataLayer(UInt size,
                                                    UInt nb_component,
                                                    const_reference value,
                                                    const ID & id)
    : ArrayBase(id) {
  allocate(size, nb_component, value);
}

/* -------------------------------------------------------------------------- */
template <typename T, ArrayAllocationType allocation_trait>
ArrayDataLayer<T, allocation_trait>::ArrayDataLayer(const ArrayDataLayer & vect,
                                                    const ID & id)
    : ArrayBase(vect, id) {
  this->data_storage = vect.data_storage;
  this->size_ = vect.size_;
  this->nb_component = vect.nb_component;
  this->values = this->data_storage.data();
}

/* -------------------------------------------------------------------------- */
template <typename T, ArrayAllocationType allocation_trait>
ArrayDataLayer<T, allocation_trait>::ArrayDataLayer(
    const std::vector<value_type> & vect) {
  this->data_storage = vect;
  this->size_ = vect.size();
  this->nb_component = 1;
  this->values = this->data_storage.data();
}

/* -------------------------------------------------------------------------- */
template <typename T, ArrayAllocationType allocation_trait>
ArrayDataLayer<T, allocation_trait> & ArrayDataLayer<T, allocation_trait>::
operator=(const ArrayDataLayer & other) {
  if (this != &other) {
    this->data_storage = other.data_storage;
    this->nb_component = other.nb_component;
    this->size_ = other.size_;
    this->values = this->data_storage.data();
  }
  return *this;
}

/* -------------------------------------------------------------------------- */
template <typename T, ArrayAllocationType allocation_trait>
ArrayDataLayer<T, allocation_trait>::ArrayDataLayer(ArrayDataLayer && other) =
    default;

/* -------------------------------------------------------------------------- */
template <typename T, ArrayAllocationType allocation_trait>
ArrayDataLayer<T, allocation_trait> & ArrayDataLayer<T, allocation_trait>::
operator=(ArrayDataLayer && other) = default;

/* -------------------------------------------------------------------------- */
template <typename T, ArrayAllocationType allocation_trait>
void ArrayDataLayer<T, allocation_trait>::allocate(UInt new_size,
                                                   UInt nb_component) {
  this->nb_component = nb_component;
  this->resize(new_size);
}

/* -------------------------------------------------------------------------- */
template <typename T, ArrayAllocationType allocation_trait>
void ArrayDataLayer<T, allocation_trait>::allocate(UInt new_size,
                                                   UInt nb_component,
                                                   const T & val) {
  this->nb_component = nb_component;
  this->resize(new_size, val);
}

/* -------------------------------------------------------------------------- */
template <typename T, ArrayAllocationType allocation_trait>
void ArrayDataLayer<T, allocation_trait>::resize(UInt new_size) {
  this->data_storage.resize(new_size * this->nb_component);
  this->values = this->data_storage.data();
  this->size_ = new_size;
}

/* -------------------------------------------------------------------------- */
template <typename T, ArrayAllocationType allocation_trait>
void ArrayDataLayer<T, allocation_trait>::resize(UInt new_size,
                                                 const T & value) {
  this->data_storage.resize(new_size * this->nb_component, value);
  this->values = this->data_storage.data();
  this->size_ = new_size;
}

/* -------------------------------------------------------------------------- */
template <typename T, ArrayAllocationType allocation_trait>
void ArrayDataLayer<T, allocation_trait>::reserve(UInt size, UInt new_size) {
  if (new_size != UInt(-1)) {
    this->data_storage.resize(new_size * this->nb_component);
  }

  this->data_storage.reserve(size * this->nb_component);
  this->values = this->data_storage.data();
}

/* -------------------------------------------------------------------------- */
/**
 * append a tuple to the array with the value value for all components
 * @param value the new last tuple or the array will contain nb_component copies
 * of value
 */
template <typename T, ArrayAllocationType allocation_trait>
inline void ArrayDataLayer<T, allocation_trait>::push_back(const T & value) {
  this->data_storage.push_back(value);
  this->values = this->data_storage.data();
  this->size_ += 1;
}

/* -------------------------------------------------------------------------- */
/**
 * append a matrix or a vector to the array
 * @param new_elem a reference to a Matrix<T> or Vector<T> */
template <typename T, ArrayAllocationType allocation_trait>
template <template <typename> class C, typename>
inline void
ArrayDataLayer<T, allocation_trait>::push_back(const C<T> & new_elem) {
  AKANTU_DEBUG_ASSERT(
      nb_component == new_elem.size(),
      "The vector("
          << new_elem.size()
          << ") as not a size compatible with the Array (nb_component="
          << nb_component << ").");
  for (UInt i = 0; i < new_elem.size(); ++i) {
    this->data_storage.push_back(new_elem[i]);
  }
  this->values = this->data_storage.data();
  this->size_ += 1;
}

/* -------------------------------------------------------------------------- */
template <typename T, ArrayAllocationType allocation_trait>
inline UInt ArrayDataLayer<T, allocation_trait>::getAllocatedSize() const {
  return this->data_storage.capacity() / this->nb_component;
}

/* -------------------------------------------------------------------------- */
template <typename T, ArrayAllocationType allocation_trait>
inline UInt ArrayDataLayer<T, allocation_trait>::getMemorySize() const {
  return this->data_storage.capacity() * sizeof(T);
}

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template <typename T>
class ArrayDataLayer<T, ArrayAllocationType::_pod> : public ArrayBase {
public:
  using value_type = T;
  using reference = value_type &;
  using pointer_type = value_type *;
  using const_reference = const value_type &;

public:
  virtual ~ArrayDataLayer() { deallocate(); }

  /// Allocation of a new vector
  ArrayDataLayer(UInt size = 0, UInt nb_component = 1, const ID & id = "")
      : ArrayBase(id) {
    allocate(size, nb_component);
  }

  /// Allocation of a new vector with a default value
  ArrayDataLayer(UInt size, UInt nb_component, const_reference value,
                 const ID & id = "")
      : ArrayBase(id) {
    allocate(size, nb_component, value);
  }

  /// Copy constructor (deep copy)
  ArrayDataLayer(const ArrayDataLayer & vect, const ID & id = "")
      : ArrayBase(vect, id) {
    allocate(vect.size(), vect.getNbComponent());
    std::copy_n(vect.storage(), this->size_ * this->nb_component, values);
  }

  /// Copy constructor (deep copy)
  explicit ArrayDataLayer(const std::vector<value_type> & vect) {
    allocate(vect.size(), 1);
    std::copy_n(vect.data(), this->size_ * this->nb_component, values);
  }

  // copy operator
  inline ArrayDataLayer & operator=(const ArrayDataLayer & other) {
    if (this != &other) {
      allocate(other.size(), other.getNbComponent());
      std::copy_n(other.storage(), this->size_ * this->nb_component, values);
    }
    return *this;
  }

  // move constructor
  inline ArrayDataLayer(ArrayDataLayer && other) = default;

  // move assign
  inline ArrayDataLayer & operator=(ArrayDataLayer && other) = default;

protected:
  // deallocate the memory
  virtual void deallocate() { free(this->values); }

  // allocate the memory
  virtual inline void allocate(UInt size, UInt nb_component) {
    if (size != 0) { // malloc can return a non NULL pointer in case size is 0
      this->values =
          static_cast<T *>(std::malloc(nb_component * size * sizeof(T)));
    }

    if (this->values == nullptr and size != 0) {
      throw std::bad_alloc();
    }
    this->nb_component = nb_component;
    this->allocated_size = this->size_ = size;
  }

  // allocate and initialize the memory
  virtual inline void allocate(UInt size, UInt nb_component, const T & value) {
    allocate(size, nb_component);
    std::fill_n(values, size * nb_component, value);
  }

public:
  /// append a tuple of size nb_component containing value
  inline void push_back(const_reference value) {
    resize(this->size_ + 1, value);
  }

  /// append a Vector or a Matrix
  template <template <typename> class C,
            typename = std::enable_if_t<aka::is_tensor<C<T>>::value>>
  inline void push_back(const C<T> & new_elem) {
    AKANTU_DEBUG_ASSERT(
        nb_component == new_elem.size(),
        "The vector("
            << new_elem.size()
            << ") as not a size compatible with the Array (nb_component="
            << nb_component << ").");
    this->resize(this->size_ + 1);
    std::copy_n(new_elem.storage(), new_elem.size(),
                values + this->nb_component * (this->size_ - 1));
  }

  /// changes the allocated size but not the size
  virtual void reserve(UInt size, UInt new_size = UInt(-1)) {
    UInt tmp_size = this->size_;
    if (new_size != UInt(-1))
      tmp_size = new_size;
    this->resize(size);
    this->size_ = std::min(this->size_, tmp_size);
  }

  /// change the size of the Array
  virtual void resize(UInt size) {
    if (size * this->nb_component == 0) {
      free(values);
      values = nullptr;
      this->allocated_size = 0;
    } else {
      if (this->values == nullptr) {
        this->allocate(size, this->nb_component);
        return;
      }

      Int diff = size - allocated_size;
      UInt size_to_allocate = (std::abs(diff) > AKANTU_MIN_ALLOCATION)
                                  ? size
                                  : (diff > 0)
                                        ? allocated_size + AKANTU_MIN_ALLOCATION
                                        : allocated_size;

      auto * tmp_ptr = reinterpret_cast<T *>(realloc(
          this->values, size_to_allocate * this->nb_component * sizeof(T)));
      if (tmp_ptr == nullptr) {
        throw std::bad_alloc();
      }

      this->values = tmp_ptr;
      this->allocated_size = size_to_allocate;
    }

    this->size_ = size;
  }

  /// change the size of the Array and initialize the values
  virtual void resize(UInt size, const T & val) {
    UInt tmp_size = this->size_;
    this->resize(size);
    if (size > tmp_size) {
      std::fill_n(values + this->nb_component * tmp_size,
                  (size - tmp_size) * this->nb_component, val);
    }
  }

  /// get the amount of space allocated in bytes
  inline UInt getMemorySize() const override final {
    return this->allocated_size * this->nb_component * sizeof(T);
  }

  /// Get the real size allocated in memory
  inline UInt getAllocatedSize() const { return this->allocated_size; }

  /// give the address of the memory allocated for this vector
  T * storage() const { return values; };

protected:
  /// allocation type agnostic  data access
  T * values{nullptr};

  UInt allocated_size{0};
};
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/* template <class T> class AllocatorMalloc : public Allocator<T> {
public:
  T * allocate(UInt size, UInt nb_component) override final {
    auto * ptr = reinterpret_cast<T *>(malloc(nb_component * size * sizeof(T)));

    if (ptr == nullptr and size != 0) {
      throw std::bad_alloc();
    }
    return ptr;
  }

  void deallocate(T * ptr, UInt size, UInt ,
                  UInt nb_component) override final {
    if (ptr) {
      if (not is_scalar<T>::value) {
        for (UInt i = 0; i < size * nb_component; ++i) {
          (ptr + i)->~T();
        }
      }
      free(ptr);
    }
  }

  std::tuple<T *, UInt> resize(UInt new_size, UInt size, UInt allocated_size,
                               UInt nb_component, T * ptr) override final {
    UInt size_to_alloc = 0;

    if (not is_scalar<T>::value and (new_size < size)) {
      for (UInt i = new_size * nb_component; i < size * nb_component; ++i) {
        (ptr + i)->~T();
      }
    }

    // free some memory
    if (new_size == 0) {
      free(ptr);
      return std::make_tuple(nullptr, 0);
    }

    if (new_size <= allocated_size) {
      if (allocated_size - new_size > AKANTU_MIN_ALLOCATION) {
        size_to_alloc = new_size;
      } else {
        return std::make_tuple(ptr, allocated_size);
      }
    } else {
      // allocate more memory
      size_to_alloc = (new_size - allocated_size < AKANTU_MIN_ALLOCATION)
                          ? allocated_size + AKANTU_MIN_ALLOCATION
                          : new_size;
    }

    auto * tmp_ptr = reinterpret_cast<T *>(
        realloc(ptr, size_to_alloc * nb_component * sizeof(T)));
    if (tmp_ptr == nullptr) {
      throw std::bad_alloc();
    }

    return std::make_tuple(tmp_ptr, size_to_alloc);
  }
};
*/

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline auto Array<T, is_scal>::operator()(UInt i, UInt j) -> reference {
  AKANTU_DEBUG_ASSERT(this->size_ > 0,
                      "The array \"" << this->id << "\" is empty");
  AKANTU_DEBUG_ASSERT((i < this->size_) && (j < this->nb_component),
                      "The value at position ["
                          << i << "," << j << "] is out of range in array \""
                          << this->id << "\"");
  return this->values[i * this->nb_component + j];
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline auto Array<T, is_scal>::operator()(UInt i, UInt j) const
    -> const_reference {
  AKANTU_DEBUG_ASSERT(this->size_ > 0,
                      "The array \"" << this->id << "\" is empty");
  AKANTU_DEBUG_ASSERT((i < this->size_) && (j < this->nb_component),
                      "The value at position ["
                          << i << "," << j << "] is out of range in array \""
                          << this->id << "\"");
  return this->values[i * this->nb_component + j];
}

template <class T, bool is_scal>
inline auto Array<T, is_scal>::operator[](UInt i) -> reference {
  AKANTU_DEBUG_ASSERT(this->size_ > 0,
                      "The array \"" << this->id << "\" is empty");
  AKANTU_DEBUG_ASSERT((i < this->size_ * this->nb_component),
                      "The value at position ["
                          << i << "] is out of range in array \"" << this->id
                          << "\"");
  return this->values[i];
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline auto Array<T, is_scal>::operator[](UInt i) const -> const_reference {
  AKANTU_DEBUG_ASSERT(this->size_ > 0,
                      "The array \"" << this->id << "\" is empty");
  AKANTU_DEBUG_ASSERT((i < this->size_ * this->nb_component),
                      "The value at position ["
                          << i << "] is out of range in array \"" << this->id
                          << "\"");
  return this->values[i];
}

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
  AKANTU_DEBUG_ASSERT((this->size_ > 0), "The array is empty");
  AKANTU_DEBUG_ASSERT((i < this->size_), "The element at position ["
                                             << i << "] is out of range (" << i
                                             << ">=" << this->size_ << ")");

  if (i != (this->size_ - 1)) {
    for (UInt j = 0; j < this->nb_component; ++j) {
      this->values[i * this->nb_component + j] =
          this->values[(this->size_ - 1) * this->nb_component + j];
    }
  }

  this->resize(this->size_ - 1);
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
  AKANTU_DEBUG_ASSERT((this->size_ == vect.size_) &&
                          (this->nb_component == vect.nb_component),
                      "The too array don't have the same sizes");

  T * a = this->values;
  T * b = vect.storage();
  for (UInt i = 0; i < this->size_ * this->nb_component; ++i) {
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
  AKANTU_DEBUG_ASSERT((this->size_ == vect.size()) &&
                          (this->nb_component == vect.nb_component),
                      "The too array don't have the same sizes");

  T * a = this->values;
  T * b = vect.storage();
  for (UInt i = 0; i < this->size_ * this->nb_component; ++i) {
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
  T * a = this->values;
  for (UInt i = 0; i < this->size_ * this->nb_component; ++i) {
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
  bool equal = this->nb_component == array.nb_component &&
               this->size_ == array.size_ && this->id == array.id;
  if (!equal)
    return false;

  if (this->values == array.storage())
    return true;
  else
    return std::equal(this->values,
                      this->values + this->size_ * this->nb_component,
                      array.storage());
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
bool Array<T, is_scal>::operator!=(const Array<T, is_scal> & array) const {
  return !operator==(array);
}

/* -------------------------------------------------------------------------- */
/**
 * set all tuples of the array to a given vector or matrix
 * @param vm Matrix or Vector to fill the array with
 */
template <class T, bool is_scal>
template <template <typename> class C, typename>
inline void Array<T, is_scal>::set(const C<T> & vm) {
  AKANTU_DEBUG_ASSERT(
      this->nb_component == vm.size(),
      "The size of the object does not match the number of components");
  for (T * it = this->values;
       it < this->values + this->nb_component * this->size_;
       it += this->nb_component) {
    std::copy_n(vm.storage(), this->nb_component, it);
  }
}
/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
void Array<T, is_scal>::append(const Array<T> & other) {
  AKANTU_DEBUG_ASSERT(
      this->nb_component == other.nb_component,
      "Cannot append an array with a different number of component");
  UInt old_size = this->size_;
  this->resize(this->size_ + other.size());

  T * tmp = this->values + this->nb_component * old_size;
  std::copy_n(other.storage(), other.size() * this->nb_component, tmp);
}

/* -------------------------------------------------------------------------- */
/* Functions Array<T, is_scal>                                                */
/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
Array<T, is_scal>::Array(UInt size, UInt nb_component, const ID & id)
    : parent(size, nb_component, id) {}

template <>
inline Array<std::string, false>::Array(UInt size, UInt nb_component,
                                        const ID & id)
    : parent(size, nb_component, "", id) {}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
Array<T, is_scal>::Array(UInt size, UInt nb_component, const_reference value,
                         const ID & id)
    : parent(size, nb_component, value, id) {}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
Array<T, is_scal>::Array(const Array & vect, const ID & id)
    : parent(vect, id) {}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
Array<T, is_scal> & Array<T, is_scal>::
operator=(const Array<T, is_scal> & other) {
  AKANTU_DEBUG_WARNING("You are copying the array "
                       << this->id << " are you sure it is on purpose");

  if (&other == this)
    return *this;

  parent::operator=(other);

  return *this;
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
Array<T, is_scal>::Array(const std::vector<T> & vect) : parent(vect) {}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal> Array<T, is_scal>::~Array() = default;

/* -------------------------------------------------------------------------- */
/**
 * search elem in the array, return  the position of the first occurrence or
 * -1 if not found
 *  @param elem the element to look for
 *  @return index of the first occurrence of elem or -1 if elem is not present
 */
template <class T, bool is_scal>
UInt Array<T, is_scal>::find(const_reference elem) const {
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
  T * it = this->values;
  UInt i = 0;
  for (; i < this->size_; ++i) {
    if (*it == elem[0]) {
      T * cit = it;
      UInt c = 0;
      for (; (c < this->nb_component) && (*cit == elem[c]); ++c, ++cit)
        ;
      if (c == this->nb_component) {
        AKANTU_DEBUG_OUT();
        return i;
      }
    }
    it += this->nb_component;
  }
  return UInt(-1);
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
template <template <typename> class C, typename>
inline UInt Array<T, is_scal>::find(const C<T> & elem) {
  AKANTU_DEBUG_ASSERT(elem.size() == this->nb_component,
                      "Cannot find an element with a wrong size ("
                          << elem.size() << ") != " << this->nb_component);
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
    if (vect.nb_component != this->nb_component)
      AKANTU_ERROR("The two arrays do not have the same number of components");

  this->resize((vect.size_ * vect.nb_component) / this->nb_component);

  std::copy_n(vect.storage(), this->size_ * this->nb_component, this->values);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <bool is_scal> class ArrayPrintHelper {
public:
  template <typename T>
  static void print_content(const Array<T> & vect, std::ostream & stream,
                            int indent) {
    std::string space(indent, AKANTU_INDENT);

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
  std::string space(indent, AKANTU_INDENT);

  std::streamsize prec = stream.precision();
  std::ios_base::fmtflags ff = stream.flags();

  stream.setf(std::ios_base::showbase);
  stream.precision(2);

  stream << space << "Array<" << debug::demangle(typeid(T).name()) << "> ["
         << std::endl;
  stream << space << " + id             : " << this->id << std::endl;
  stream << space << " + size           : " << this->size_ << std::endl;
  stream << space << " + nb_component   : " << this->nb_component << std::endl;
  stream << space << " + allocated size : " << this->getAllocatedSize()
         << std::endl;
  stream << space
         << " + memory size    : " << printMemorySize<T>(this->getMemorySize())
         << std::endl;
  if (!AKANTU_DEBUG_LEVEL_IS_TEST())
    stream << space << " + address        : " << std::hex << this->values
           << std::dec << std::endl;

  stream.precision(prec);
  stream.flags(ff);

  if (AKANTU_DEBUG_TEST(dblDump) || AKANTU_DEBUG_LEVEL_IS_TEST()) {
    ArrayPrintHelper<is_scal or std::is_enum<T>::value>::print_content(
        *this, stream, indent);
  }

  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
/* Inline Functions ArrayBase                                                */
/* -------------------------------------------------------------------------- */

inline void ArrayBase::empty() { this->size_ = 0; }

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
  using const_reference = const R &;
  using internal_value_type = IR;
  using internal_pointer = IR *;
  using difference_type = std::ptrdiff_t;
  using iterator_category = std::random_access_iterator_tag;
  static_assert(not is_tensor, "Cannot handle tensors");

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
/**
 * Specialization for scalar types
 */
template <class T, bool is_scal>
template <class R, class daughter, class IR>
class Array<T, is_scal>::iterator_internal<R, daughter, IR, true> {
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
  using pointer_type = typename Array<T, is_scal>::pointer_type;

public:
  iterator_internal() = default;

  iterator_internal(pointer_type data, UInt _offset)
      : _offset(_offset), initial(data), ret(nullptr), ret_ptr(data) {
    AKANTU_ERROR(
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
/* Iterators                                                                  */
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

  template <typename P,
            typename = std::enable_if_t<not aka::is_tensor<P>::value>>
  const_iterator(P * data) : parent(data) {}

  template <typename UP_P, typename = std::enable_if_t<aka::is_tensor<
                               typename UP_P::element_type>::value>>
  const_iterator(UP_P && tensor) : parent(std::forward<UP_P>(tensor)) {}

  const_iterator & operator=(const const_iterator & it) = default;
};

/* -------------------------------------------------------------------------- */
template <class T, class R, bool is_tensor_ = aka::is_tensor<R>::value>
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

/* -------------------------------------------------------------------------- */
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

  template <typename P,
            typename = std::enable_if_t<not aka::is_tensor<P>::value>>
  iterator(P * data) : parent(data) {}

  template <typename UP_P, typename = std::enable_if_t<aka::is_tensor<
                               typename UP_P::element_type>::value>>
  iterator(UP_P && tensor) : parent(std::forward<UP_P>(tensor)) {}

  iterator & operator=(const iterator & it) = default;

  operator const_iterator<R>() {
    return ConstConverterIteratorHelper<T, R>::convert(*this);
  }
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

  template <typename... V> constexpr auto product_all(V &&... v) {
    std::common_type_t<int, V...> result = 1;
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
          "The iterator on "
              << debug::demangle(typeid(Arr).name())
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
  return detail::get_iterator(*this, this->values, std::forward<Ns>(ns)...,
                              this->size_);
}

template <class T, bool is_scal>
template <typename... Ns>
inline decltype(auto) Array<T, is_scal>::end(Ns &&... ns) {
  return detail::get_iterator(*this,
                              this->values + this->nb_component * this->size_,
                              std::forward<Ns>(ns)..., this->size_);
}

template <class T, bool is_scal>
template <typename... Ns>
inline decltype(auto) Array<T, is_scal>::begin(Ns &&... ns) const {
  return detail::get_iterator(*this, this->values, std::forward<Ns>(ns)...,
                              this->size_);
}

template <class T, bool is_scal>
template <typename... Ns>
inline decltype(auto) Array<T, is_scal>::end(Ns &&... ns) const {
  return detail::get_iterator(*this,
                              this->values + this->nb_component * this->size_,
                              std::forward<Ns>(ns)..., this->size_);
}

template <class T, bool is_scal>
template <typename... Ns>
inline decltype(auto) Array<T, is_scal>::begin_reinterpret(Ns &&... ns) {
  return detail::get_iterator(*this, this->values, std::forward<Ns>(ns)...);
}

template <class T, bool is_scal>
template <typename... Ns>
inline decltype(auto) Array<T, is_scal>::end_reinterpret(Ns &&... ns) {
  return detail::get_iterator(
      *this, this->values + detail::product_all(std::forward<Ns>(ns)...),
      std::forward<Ns>(ns)...);
}

template <class T, bool is_scal>
template <typename... Ns>
inline decltype(auto) Array<T, is_scal>::begin_reinterpret(Ns &&... ns) const {
  return detail::get_iterator(*this, this->values, std::forward<Ns>(ns)...);
}

template <class T, bool is_scal>
template <typename... Ns>
inline decltype(auto) Array<T, is_scal>::end_reinterpret(Ns &&... ns) const {
  return detail::get_iterator(
      *this, this->values + detail::product_all(std::forward<Ns>(ns)...),
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
        : array(array), sizes(std::move(ns)...) {}

    ArrayView(ArrayView && array_view) = default;

    ArrayView & operator=(const ArrayView & array_view) = default;
    ArrayView & operator=(ArrayView && array_view) = default;

    decltype(auto) begin() {
      return aka::apply(
          [&](auto &&... ns) { return array.get().begin_reinterpret(ns...); },
          sizes);
    }

    decltype(auto) begin() const {
      return aka::apply(
          [&](auto &&... ns) { return array.get().begin_reinterpret(ns...); },
          sizes);
    }

    decltype(auto) end() {
      return aka::apply(
          [&](auto &&... ns) { return array.get().end_reinterpret(ns...); },
          sizes);
    }

    decltype(auto) end() const {
      return aka::apply(
          [&](auto &&... ns) { return array.get().end_reinterpret(ns...); },
          sizes);
    }

    decltype(auto) size() const {
      return std::get<std::tuple_size<tuple>::value - 1>(sizes);
    }

    decltype(auto) dims() const { return std::tuple_size<tuple>::value - 1; }

  private:
    std::reference_wrapper<std::remove_reference_t<Array>> array;
    tuple sizes;
  };
} // namespace detail

/* -------------------------------------------------------------------------- */
template <typename Array, typename... Ns>
decltype(auto) make_view(Array && array, Ns... ns) {
  static_assert(aka::conjunction<std::is_integral<std::decay_t<Ns>>...>::value,
                "Ns should be integral types");
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
inline typename Array<T, is_scal>::template iterator<R>
Array<T, is_scal>::erase(const iterator<R> & it) {
  T * curr = it.data();
  UInt pos = (curr - this->values) / this->nb_component;
  erase(pos);
  iterator<R> rit = it;
  return --rit;
}

} // namespace akantu

#endif /* __AKANTU_AKA_ARRAY_TMPL_HH__ */
