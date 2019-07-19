/**
 * @file   aka_array.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Jan 16 2018
 *
 * @brief  Array container for Akantu
 * This container differs from the std::vector from the fact it as 2 dimensions
 * a main dimension and the size stored per entries
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

#ifndef __AKANTU_VECTOR_HH__
#define __AKANTU_VECTOR_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */
#include <typeinfo>
#include <vector>
/* -------------------------------------------------------------------------- */

namespace akantu {

/// class that afford to store vectors in static memory
class ArrayBase {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  explicit ArrayBase(ID id = "") : id(std::move(id)) {}
  ArrayBase(const ArrayBase & other, const ID & id = "") {
    this->id = (id == "") ? other.id : id;
  }

  ArrayBase(ArrayBase && other) = default;
  ArrayBase & operator=(const ArrayBase & other) = default;
  // ArrayBase & operator=(ArrayBase && other) = default;

  virtual ~ArrayBase() = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// get the amount of space allocated in bytes
  virtual UInt getMemorySize() const = 0;

  /// set the size to zero without freeing the allocated space
  inline void empty();

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const = 0;

  /* ------------------------------------------------------------------------ */
  /* Accessors */
  /* ------------------------------------------------------------------------ */
public:
  /// Get the Size of the Array
  UInt size() const { return size_; }
  /// Get the number of components
  AKANTU_GET_MACRO(NbComponent, nb_component, UInt);
  /// Get the name of th array
  AKANTU_GET_MACRO(ID, id, const ID &);
  /// Set the name of th array
  AKANTU_SET_MACRO(ID, id, const ID &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// id of the vector
  ID id{""};

  /// the size used
  UInt size_{0};

  /// number of components
  UInt nb_component{1};
};

/* -------------------------------------------------------------------------- */
namespace {
  template <std::size_t dim, typename T> struct IteratorHelper {};

  template <typename T> struct IteratorHelper<0, T> { using type = T; };
  template <typename T> struct IteratorHelper<1, T> { using type = Vector<T>; };
  template <typename T> struct IteratorHelper<2, T> { using type = Matrix<T>; };
  template <typename T> struct IteratorHelper<3, T> {
    using type = Tensor3<T>;
  };

  template <std::size_t dim, typename T>
  using IteratorHelper_t = typename IteratorHelper<dim, T>::type;
} // namespace

/* -------------------------------------------------------------------------- */
/* Memory handling layer                                                      */
/* -------------------------------------------------------------------------- */
enum class ArrayAllocationType {
  _default,
  _pod,
};

template <typename T>
struct ArrayAllocationTrait
    : public std::conditional_t<
          std::is_scalar<T>::value,
          std::integral_constant<ArrayAllocationType,
                                 ArrayAllocationType::_pod>,
          std::integral_constant<ArrayAllocationType,
                                 ArrayAllocationType::_default>> {};

/* -------------------------------------------------------------------------- */
template <typename T,
          ArrayAllocationType allocation_trait = ArrayAllocationTrait<T>::value>
class ArrayDataLayer : public ArrayBase {
public:
  using value_type = T;
  using reference = value_type &;
  using pointer_type = value_type *;
  using const_reference = const value_type &;

public:
  virtual ~ArrayDataLayer() = default;

  /// Allocation of a new vector
  explicit ArrayDataLayer(UInt size = 0, UInt nb_component = 1,
                          const ID & id = "");

  /// Allocation of a new vector with a default value
  ArrayDataLayer(UInt size, UInt nb_component, const_reference value,
                 const ID & id = "");

  /// Copy constructor (deep copy)
  ArrayDataLayer(const ArrayDataLayer & vect, const ID & id = "");

  /// Copy constructor (deep copy)
  explicit ArrayDataLayer(const std::vector<value_type> & vect);

  // copy operator
  ArrayDataLayer & operator=(const ArrayDataLayer & other);

  // move constructor
  ArrayDataLayer(ArrayDataLayer && other);

  // move assign
  ArrayDataLayer & operator=(ArrayDataLayer && other);

protected:
  // deallocate the memory
  virtual void deallocate() {}

  // allocate the memory
  virtual void allocate(UInt size, UInt nb_component);

  // allocate and initialize the memory
  virtual void allocate(UInt size, UInt nb_component, const T & value);

public:
  /// append a tuple of size nb_component containing value
  inline void push_back(const_reference value);
  /// append a vector
  // inline void push_back(const value_type new_elem[]);

  /// append a Vector or a Matrix
  template <template <typename> class C,
            typename = std::enable_if_t<aka::is_tensor<C<T>>::value>>
  inline void push_back(const C<T> & new_elem);

  /// changes the allocated size but not the size, if new_size = 0, the size is
  /// set to min(current_size and reserve size)
  virtual void reserve(UInt size, UInt new_size = UInt(-1));

  /// change the size of the Array
  virtual void resize(UInt size);

  /// change the size of the Array and initialize the values
  virtual void resize(UInt size, const T & val);

  /// get the amount of space allocated in bytes
  inline UInt getMemorySize() const override;

  /// Get the real size allocated in memory
  inline UInt getAllocatedSize() const;

  /// give the address of the memory allocated for this vector
  T * storage() const { return values; };

protected:
  /// allocation type agnostic  data access
  T * values{nullptr};

  /// data storage
  std::vector<T> data_storage;
};

/* -------------------------------------------------------------------------- */
/* Actual Array                                                               */
/* -------------------------------------------------------------------------- */
template <typename T, bool is_scal> class Array : public ArrayDataLayer<T> {
private:
  using parent = ArrayDataLayer<T>;
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  using value_type = typename parent::value_type;
  using reference = typename parent::reference;
  using pointer_type = typename parent::pointer_type;
  using const_reference = typename parent::const_reference;

  ~Array() override;

  Array() : Array(0){};

  /// Allocation of a new vector
  explicit Array(UInt size, UInt nb_component = 1, const ID & id = "");

  /// Allocation of a new vector with a default value
  explicit Array(UInt size, UInt nb_component, const_reference value,
                 const ID & id = "");

  /// Copy constructor (deep copy if deep=true)
  Array(const Array & vect, const ID & id = "");

  /// Copy constructor (deep copy)
  explicit Array(const std::vector<T> & vect);

  // copy operator
  Array & operator=(const Array & other);

  // move constructor
  Array(Array && other) = default;

  // move assign
  Array & operator=(Array && other) = default;

  /* ------------------------------------------------------------------------ */
  /* Iterator                                                                 */
  /* ------------------------------------------------------------------------ */
  /// \todo protected: does not compile with intel  check why
public:
  template <class R, class it, class IR = R,
            bool is_tensor_ = aka::is_tensor<std::decay_t<R>>::value>
  class iterator_internal;

public:
  /* ------------------------------------------------------------------------ */

  /* ------------------------------------------------------------------------ */
  template <typename R = T> class const_iterator;
  template <typename R = T> class iterator;

  /* ------------------------------------------------------------------------ */

  /// iterator for Array of nb_component = 1
  using scalar_iterator = iterator<T>;
  /// const_iterator for Array of nb_component = 1
  using const_scalar_iterator = const_iterator<T>;

  /// iterator returning Vectors of size n  on entries of Array with
  /// nb_component = n
  using vector_iterator = iterator<Vector<T>>;
  /// const_iterator returning Vectors of n size on entries of Array with
  /// nb_component = n
  using const_vector_iterator = const_iterator<Vector<T>>;

  /// iterator returning Matrices of size (m, n) on entries of Array with
  /// nb_component = m*n
  using matrix_iterator = iterator<Matrix<T>>;
  /// const iterator returning Matrices of size (m, n) on entries of Array with
  /// nb_component = m*n
  using const_matrix_iterator = const_iterator<Matrix<T>>;

  /// iterator returning Tensor3 of size (m, n, k) on entries of Array with
  /// nb_component = m*n*k
  using tensor3_iterator = iterator<Tensor3<T>>;
  /// const iterator returning Tensor3 of size (m, n, k) on entries of Array
  /// with nb_component = m*n*k
  using const_tensor3_iterator = const_iterator<Tensor3<T>>;

  /* ------------------------------------------------------------------------ */
  template <typename... Ns> inline decltype(auto) begin(Ns &&... n);
  template <typename... Ns> inline decltype(auto) end(Ns &&... n);

  template <typename... Ns> inline decltype(auto) begin(Ns &&... n) const;
  template <typename... Ns> inline decltype(auto) end(Ns &&... n) const;

  template <typename... Ns> inline decltype(auto) begin_reinterpret(Ns &&... n);
  template <typename... Ns> inline decltype(auto) end_reinterpret(Ns &&... n);

  template <typename... Ns>
  inline decltype(auto) begin_reinterpret(Ns &&... n) const;
  template <typename... Ns>
  inline decltype(auto) end_reinterpret(Ns &&... n) const;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// search elem in the vector, return  the position of the first occurrence or
  /// -1 if not found
  UInt find(const_reference elem) const;

  /// @see Array::find(const_reference elem) const
  UInt find(T elem[]) const;

  inline void push_back(const_reference value) { parent::push_back(value); }

  /// append a Vector or a Matrix
  template <template <typename> class C,
            typename = std::enable_if_t<aka::is_tensor<C<T>>::value>>
  inline void push_back(const C<T> & new_elem) {
    parent::push_back(new_elem);
  }

  template <typename Ret> inline void push_back(const iterator<Ret> & it) {
    push_back(*it);
  }

  /// erase the value at position i
  inline void erase(UInt i);
  /// ask Nico, clarify
  template <typename R> inline iterator<R> erase(const iterator<R> & it);

  /// @see Array::find(const_reference elem) const
  template <template <typename> class C,
            typename = std::enable_if_t<aka::is_tensor<C<T>>::value>>
  inline UInt find(const C<T> & elem);

  /// set all entries of the array to the value t
  /// @param t value to fill the array with
  inline void set(T t) {
    std::fill_n(this->values, this->size_ * this->nb_component, t);
  }

  /// set all entries of the array to 0
  inline void clear() { set(T()); }

  /// set all tuples of the array to a given vector or matrix
  /// @param vm Matrix or Vector to fill the array with
  template <template <typename> class C,
            typename = std::enable_if_t<aka::is_tensor<C<T>>::value>>
  inline void set(const C<T> & vm);

  /// Append the content of the other array to the current one
  void append(const Array<T> & other);

  /// copy another Array in the current Array, the no_sanity_check allows you to
  /// force the copy in cases where you know what you do with two non matching
  /// Arrays in terms of n
  void copy(const Array<T, is_scal> & other, bool no_sanity_check = false);

  /// function to print the containt of the class
  void printself(std::ostream & stream, int indent = 0) const override;

  /* ------------------------------------------------------------------------ */
  /* Operators                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// substraction entry-wise
  Array<T, is_scal> & operator-=(const Array<T, is_scal> & other);
  /// addition entry-wise
  Array<T, is_scal> & operator+=(const Array<T, is_scal> & other);
  /// multiply evry entry by alpha
  Array<T, is_scal> & operator*=(const T & alpha);

  /// check if the array are identical entry-wise
  bool operator==(const Array<T, is_scal> & other) const;
  /// @see Array::operator==(const Array<T, is_scal> & other) const
  bool operator!=(const Array<T, is_scal> & other) const;

  /// return a reference to the j-th entry of the i-th tuple
  inline reference operator()(UInt i, UInt j = 0);
  /// return a const reference to the j-th entry of the i-th tuple
  inline const_reference operator()(UInt i, UInt j = 0) const;

  /// return a reference to the ith component of the 1D array
  inline reference operator[](UInt i);
  /// return a const reference to the ith component of the 1D array
  inline const_reference operator[](UInt i) const;
};

/* -------------------------------------------------------------------------- */
/* Inline Functions Array<T, is_scal>                                         */
/* -------------------------------------------------------------------------- */
template <typename T, bool is_scal>
inline std::ostream & operator<<(std::ostream & stream,
                                 const Array<T, is_scal> & _this) {
  _this.printself(stream);
  return stream;
}

/* -------------------------------------------------------------------------- */
/* Inline Functions ArrayBase                                                 */
/* -------------------------------------------------------------------------- */
inline std::ostream & operator<<(std::ostream & stream,
                                 const ArrayBase & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#include "aka_array_tmpl.hh"
#include "aka_types.hh"

#endif /* __AKANTU_VECTOR_HH__ */
