/**
 * @file   aka_array.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Fri Jan 22 2016
 *
 * @brief  Array container for Akantu
 * This container differs from the std::vector from the fact it as 2 dimensions
 * a main dimension and the size stored per entries
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

#ifndef __AKANTU_VECTOR_HH__
#define __AKANTU_VECTOR_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */
#include <typeinfo>
#include <vector>
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/// class that afford to store vectors in static memory
class ArrayBase {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  ArrayBase(const ID & id = "");

  virtual ~ArrayBase();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// get the amount of space allocated in bytes
  inline UInt getMemorySize() const;

  /// set the size to zero without freeing the allocated space
  inline void empty();

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors */
  /* ------------------------------------------------------------------------ */
public:
  /// Get the real size allocated in memory
  AKANTU_GET_MACRO(AllocatedSize, allocated_size, UInt);
  /// Get the Size of the Array
  AKANTU_GET_MACRO(Size, size, UInt);
  /// Get the number of components
  AKANTU_GET_MACRO(NbComponent, nb_component, UInt);
  /// Get the name of th array
  AKANTU_GET_MACRO(ID, id, const ID &);
  /// Set the name of th array
  AKANTU_SET_MACRO(ID, id, const ID &);

  // AKANTU_GET_MACRO(Tag, tag, const std::string &);
  // AKANTU_SET_MACRO(Tag, tag, const std::string &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// id of the vector
  ID id;

  /// the size allocated
  UInt allocated_size;

  /// the size used
  UInt size;

  /// number of components
  UInt nb_component;

  /// size of the stored type
  UInt size_of_type;

  // /// User defined tag
  // std::string tag;
};

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template <typename T, bool is_scal> class Array : public ArrayBase {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  typedef T value_type;
  typedef value_type & reference;
  typedef value_type * pointer_type;
  typedef const value_type & const_reference;

  /// Allocation of a new vector
  inline Array(UInt size = 0, UInt nb_component = 1, const ID & id = "");

  /// Allocation of a new vector with a default value
  Array(UInt size, UInt nb_component, const value_type def_values[],
        const ID & id = "");

  /// Allocation of a new vector with a default value
  Array(UInt size, UInt nb_component, const_reference value,
        const ID & id = "");

  /// Copy constructor (deep copy if deep=true)
  Array(const Array<value_type, is_scal> & vect, bool deep = true,
        const ID & id = "");

#ifndef SWIG
  /// Copy constructor (deep copy)
  Array(const std::vector<value_type> & vect);
#endif

  virtual inline ~Array();

  Array & operator=(const Array & a) {
    /// this is to let STL allocate and copy arrays in the case of
    /// std::vector::resize
    AKANTU_DEBUG_ASSERT(this->size == 0, "Cannot copy akantu::Array");
    return const_cast<Array &>(a);
  }

  /* ------------------------------------------------------------------------ */
  /* Iterator                                                                 */
  /* ------------------------------------------------------------------------ */
  /// \todo protected: does not compile with intel  check why
public:
  template <class R, class IR = R, bool issame = is_same<IR, T>::value>
  class iterator_internal;

public:
  /* ------------------------------------------------------------------------ */

  /* ------------------------------------------------------------------------ */
  template <typename R = T> class const_iterator;
  template <typename R = T> class iterator;

  /* ------------------------------------------------------------------------ */

  /// iterator for Array of nb_component = 1
  typedef iterator<T> scalar_iterator;
  /// const_iterator for Array of nb_component = 1
  typedef const_iterator<T> const_scalar_iterator;

  /// iterator rerturning Vectors of size n  on entries of Array with
  /// nb_component = n
  typedef iterator<Vector<T> > vector_iterator;
  /// const_iterator rerturning Vectors of n size on entries of Array with
  /// nb_component = n
  typedef const_iterator<Vector<T> > const_vector_iterator;

  /// iterator rerturning Matrices of size (m, n) on entries of Array with
  /// nb_component = m*n
  typedef iterator<Matrix<T> > matrix_iterator;
  /// const iterator rerturning Matrices of size (m, n) on entries of Array with
  /// nb_component = m*n
  typedef const_iterator<Matrix<T> > const_matrix_iterator;

  /* ------------------------------------------------------------------------ */

  /// Get an iterator that behaves like a pointer T * to the first entry
  inline scalar_iterator begin();
  /// Get an iterator that behaves like a pointer T * to the end of the Array
  inline scalar_iterator end();
  /// Get a const_iterator to the beginging of an Array of scalar
  inline const_scalar_iterator begin() const;
  /// Get a const_iterator to the end of an Array of scalar
  inline const_scalar_iterator end() const;

  /// Get a scalar_iterator on the beginning of the Array considered of shape (new_size)
  inline scalar_iterator begin_reinterpret(UInt new_size);
  /// Get a scalar_iterator on the end of the Array considered of shape (new_size)
  inline scalar_iterator end_reinterpret(UInt new_size);
  /// Get a const_scalar_iterator on the beginning of the Array considered of shape (new_size)
  inline const_scalar_iterator begin_reinterpret(UInt new_size) const;
  /// Get a const_scalar_iterator on the end of the Array considered of shape (new_size)
  inline const_scalar_iterator end_reinterpret(UInt new_size) const;

  /* ------------------------------------------------------------------------ */
  /// Get a vector_iterator on the beginning of the Array
  inline vector_iterator begin(UInt n);
  /// Get a vector_iterator on the end of the Array
  inline vector_iterator end(UInt n);
  /// Get a vector_iterator on the beginning of the Array
  inline const_vector_iterator begin(UInt n) const;
  /// Get a vector_iterator on the end of the Array
  inline const_vector_iterator end(UInt n) const;

  /// Get a vector_iterator on the begining of the Array considered of shape
  /// (new_size, n)
  inline vector_iterator begin_reinterpret(UInt n, UInt new_size);
  /// Get a vector_iterator on the end of the Array considered of shape
  /// (new_size, n)
  inline vector_iterator end_reinterpret(UInt n, UInt new_size);
  /// Get a const_vector_iterator on the begining of the Array considered of
  /// shape (new_size, n)
  inline const_vector_iterator begin_reinterpret(UInt n, UInt new_size) const;
  /// Get a const_vector_iterator on the end of the Array considered of shape
  /// (new_size, n)
  inline const_vector_iterator end_reinterpret(UInt n, UInt new_size) const;

  /* ------------------------------------------------------------------------ */
  /// Get a matrix_iterator on the begining of the Array (Matrices of size (m,
  /// n))
  inline matrix_iterator begin(UInt m, UInt n);
  /// Get a matrix_iterator on the end of the Array (Matrices of size (m, n))
  inline matrix_iterator end(UInt m, UInt n);
  /// Get a const_matrix_iterator on the begining of the Array (Matrices of size
  /// (m, n))
  inline const_matrix_iterator begin(UInt m, UInt n) const;
  /// Get a const_matrix_iterator on the end of the Array (Matrices of size (m,
  /// n))
  inline const_matrix_iterator end(UInt m, UInt n) const;

  /// Get a matrix_iterator on the begining of the Array considered of shape
  /// (new_size, m*n)
  inline matrix_iterator begin_reinterpret(UInt m, UInt n, UInt size);
  /// Get a matrix_iterator on the end of the Array considered of shape
  /// (new_size, m*n)
  inline matrix_iterator end_reinterpret(UInt m, UInt n, UInt size);
  /// Get a const_matrix_iterator on the begining of the Array considered of
  /// shape (new_size, m*n)
  inline const_matrix_iterator begin_reinterpret(UInt m, UInt n,
                                                 UInt size) const;
  /// Get a const_matrix_iterator on the end of the Array considered of shape
  /// (new_size, m*n)
  inline const_matrix_iterator end_reinterpret(UInt m, UInt n, UInt size) const;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// append a tuple of size nb_component containing value
  inline void push_back(const_reference value);
  /// append a vector
  //inline void push_back(const value_type new_elem[]);

  /// append a Vector or a Matrix
  template <template <typename> class C>
  inline void push_back(const C<T> & new_elem);
  /// append the value of the iterator
  template <typename Ret> inline void push_back(const iterator<Ret> & it);

  /// erase the value at position i
  inline void erase(UInt i);
  /// ask Nico, clarify
  template <typename R> inline iterator<R> erase(const iterator<R> & it);

  /// change the size of the Array
  virtual void resize(UInt size);

  /// change the size of the Array and initialize the values
  virtual void resize(UInt size, const T & val);

  /// change the number of components by interlacing data
  /// @param multiplicator number of interlaced components add
  /// @param block_size blocks of data in the array
  /// Examaple for block_size = 2, multiplicator = 2
  /// array = oo oo oo -> new array = oo nn nn oo nn nn oo nn nn
  void extendComponentsInterlaced(UInt multiplicator, UInt stride);

  /// search elem in the vector, return  the position of the first occurrence or
  /// -1 if not found
  Int find(const_reference elem)
      const; /// @see Array::find(const_reference elem) const
  Int find(T elem[]) const;
  /// @see Array::find(const_reference elem) const
  template <template <typename> class C> inline Int find(const C<T> & elem);

  /// set all entries of the array to 0
  inline void clear() { std::fill_n(values, size * nb_component, T()); }

  /// set all entries of the array to the value t
  /// @param t value to fill the array with
  inline void set(T t) { std::fill_n(values, size * nb_component, t); }

  /// set all tuples of the array to a given vector or matrix
  /// @param vm Matrix or Vector to fill the array with
  template <template <typename> class C> inline void set(const C<T> & vm);

  /// Append the content of the other array to the current one
  void append(const Array<T> & other);

  /// copy another Array in the current Array, the no_sanity_check allows you to
  /// force the copy in cases where you know what you do with two non matching
  /// Arrays in terms of n
  void copy(const Array<T, is_scal> & other, bool no_sanity_check = false);

  /// give the address of the memory allocated for this vector
  T * storage() const { return values; };

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

protected:
  /// perform the allocation for the constructors
  void allocate(UInt size, UInt nb_component = 1);

  /// resize initializing with uninitialized_fill if fill is set
  void resizeUnitialized(UInt new_size, bool fill, const T & val = T());

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

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// array of values
  T * values; // /!\ very dangerous
};

#include "aka_array_tmpl.hh"

__END_AKANTU__

#include "aka_types.hh"

__BEGIN_AKANTU__

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

__END_AKANTU__

#endif /* __AKANTU_VECTOR_HH__ */
