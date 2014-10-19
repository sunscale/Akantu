/**
 * @file   aka_array_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Jul 15 2010
 * @date last modification: Tue Sep 02 2014
 *
 * @brief  Inline functions of the classes Array<T> and ArrayBase
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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

__END_AKANTU__

#include <memory>

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline T & Array<T, is_scal>::operator()(UInt i, UInt j) {
  AKANTU_DEBUG_ASSERT(size > 0,
		      "The array \"" << id << "\" is empty");
  AKANTU_DEBUG_ASSERT((i < size) && (j < nb_component),
		      "The value at position [" << i << "," << j
		      << "] is out of range in array \"" << id << "\"");
  return values[i*nb_component + j];
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline const T & Array<T, is_scal>::operator()(UInt i, UInt j) const {
  AKANTU_DEBUG_ASSERT(size > 0,
		      "The array \"" << id << "\" is empty");
  AKANTU_DEBUG_ASSERT((i < size) && (j < nb_component),
		      "The value at position [" << i << "," << j
		      << "] is out of range in array \"" << id << "\"");
  return values[i*nb_component + j];
}

template <class T, bool is_scal>
inline T & Array<T, is_scal>::operator[](UInt i) {
  AKANTU_DEBUG_ASSERT(size > 0,
                      "The array \"" << id << "\" is empty");
  AKANTU_DEBUG_ASSERT((i < size*nb_component),
                      "The value at position [" << i << "] is out of range in array \"" << id << "\"");
  return values[i];
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline const T & Array<T, is_scal>::operator[](UInt i) const {
  AKANTU_DEBUG_ASSERT(size > 0,
                      "The array \"" << id << "\" is empty");
  AKANTU_DEBUG_ASSERT((i < size*nb_component),
                      "The value at position [" << i << "] is out of range in array \"" << id << "\"");
  return values[i];
}



/* -------------------------------------------------------------------------- */
/**
 * append a tuple to the array with the value value for all
 * components
 * @param value the new last tuple or the array will contain nb_component copies of value
 */
template <class T, bool is_scal>
inline void Array<T, is_scal>::push_back(const T & value) {
  UInt pos = size;

  resizeUnitialized(size+1);

  std::uninitialized_fill_n(values + pos * nb_component, nb_component, value);
}

/* -------------------------------------------------------------------------- */
/**
 * append a tuple to the array
 * @param new_elem a C-array containing the values to be copied to the end of the array */
template <class T, bool is_scal>
inline void Array<T, is_scal>::push_back(const T new_elem[]) {
  UInt pos = size;

  resizeUnitialized(size+1);

  T * tmp = values + nb_component * pos;
  std::uninitialized_copy(new_elem, new_elem + nb_component, tmp);
}

/* -------------------------------------------------------------------------- */
/**
 * append a matrix or a vector to the array
 * @param new_elem a reference to a Matrix<T> or Vector<T> */
template <class T, bool is_scal>
template<template<typename> class C>
inline void Array<T, is_scal>::push_back(const C<T> & new_elem) {
  AKANTU_DEBUG_ASSERT(nb_component == new_elem.size(),
		      "The vector("<< new_elem.size() <<") as not a size compatible with the Array (nb_component=" << nb_component << ").");
  UInt pos = size;
  resizeUnitialized(size+1);

  T * tmp = values + nb_component * pos;
  std::uninitialized_copy(new_elem.storage(), new_elem.storage() + nb_component, tmp);
}

/* -------------------------------------------------------------------------- */
/**
 * append a tuple to the array
 * @param it an iterator to the tuple to be copied to the end of the array */
template <class T, bool is_scal>
template<class Ret>
inline void Array<T, is_scal>::push_back(const Array<T, is_scal>::iterator<Ret> & it) {
  UInt pos = size;

  resizeUnitialized(size+1);

  T * tmp = values + nb_component * pos;
  T * new_elem = it.data();
  std::uninitialized_copy(new_elem, new_elem + nb_component, tmp);
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
template <class T, bool is_scal>
inline void Array<T, is_scal>::erase(UInt i){
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT((size > 0),
		      "The array is empty");
  AKANTU_DEBUG_ASSERT((i < size),
		      "The element at position [" << i << "] is out of range (" << i << ">=" << size << ")");


  if(i != (size - 1)) {
    for (UInt j = 0; j < nb_component; ++j) {
      values[i*nb_component + j] = values[(size-1)*nb_component + j];
    }
  }

  resize(size - 1);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Subtract another array entry by entry from this array in place. Both arrays must
 * have the same size and nb_component. If the arrays have different shapes,
 * code compiled in debug mode will throw an expeption and optimised code
 * will behave in an unpredicted manner
 * @param other array to subtract from this
 * @return reference to modified this
 */
template <class T, bool is_scal>
Array<T, is_scal> & Array<T, is_scal>::operator-=(const Array<T, is_scal> & vect) {
  AKANTU_DEBUG_ASSERT((size == vect.size) && (nb_component == vect.nb_component),
		      "The too array don't have the same sizes");

  T * a = values;
  T * b = vect.storage();
  for (UInt i = 0; i < size*nb_component; ++i) {
    *a -= *b;
    ++a;++b;
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
Array<T, is_scal> & Array<T, is_scal>::operator+=(const Array<T, is_scal> & vect) {
  AKANTU_DEBUG_ASSERT((size == vect.size) && (nb_component == vect.nb_component),
		      "The too array don't have the same sizes");

  T * a = values;
  T * b = vect.storage();
  for (UInt i = 0; i < size*nb_component; ++i) {
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
  for (UInt i = 0; i < size*nb_component; ++i) {
    *a++ *= alpha;
  }

  return *this;
}

/* -------------------------------------------------------------------------- */
/**
 * Compare this array element by element to another.
 * @param other array to compare to
 * @return true it all element are equal and arrays have the same shape, else false
 */
template <class T, bool is_scal>
bool Array<T, is_scal>::operator==(const Array<T, is_scal> & array) const {
  bool equal = nb_component == array.nb_component && size == array.size && id == array.id;
  if(!equal) return false;

  if(values == array.storage()) return true;
  else return std::equal(values, values + size*nb_component,
                         array.storage());
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
bool Array<T, is_scal>::operator!=(const Array<T, is_scal> & array) const {
  return !operator==(array);
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
template<template<typename> class C>
inline void Array<T, is_scal>::set(const C<T> & vm) {
  AKANTU_DEBUG_ASSERT(nb_component == vm.size(),
		      "The size of the object does not match the number of components");
  for (T * it = values;
       it < values + nb_component * size;
       it += nb_component) {
    std::copy(vm.storage(), vm.storage() + nb_component, it);
  }
}

/* -------------------------------------------------------------------------- */
/* Functions Array<T, is_scal>                                               */
/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
Array<T, is_scal>::Array (UInt size,
                          UInt nb_component,
                          const ID & id) :
  ArrayBase(id), values(NULL) {
  AKANTU_DEBUG_IN();
  allocate(size, nb_component);

  if(!is_scal) {
    T val = T();
    std::uninitialized_fill(values, values + size*nb_component, val);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
Array<T, is_scal>::Array (UInt size,
                          UInt nb_component,
                          const T def_values[],
                          const ID & id) :
  ArrayBase(id), values(NULL) {
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
Array<T, is_scal>::Array (UInt size,
                          UInt nb_component,
                          const T & value,
                          const ID & id) :
  ArrayBase(id), values(NULL) {
  AKANTU_DEBUG_IN();
  allocate(size, nb_component);

  std::uninitialized_fill_n(values, size*nb_component, value);

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
Array<T, is_scal>::Array(const Array<T, is_scal> & vect,
                         bool deep,
                         const ID & id) {
  AKANTU_DEBUG_IN();
  this->id = (id == "") ? vect.id : id;

  if (deep) {
    allocate(vect.size, vect.nb_component);
    T * tmp = values;
    std::uninitialized_copy(vect.storage(), vect.storage() + size * nb_component, tmp);
  } else {
    this->values = vect.storage();
    this->size = vect.size;
    this->nb_component = vect.nb_component;
    this->allocated_size = vect.allocated_size;
    this->size_of_type = vect.size_of_type;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
Array<T, is_scal>::Array(const std::vector<T>& vect) {
  AKANTU_DEBUG_IN();
  this->id = "";

  allocate(vect.size(), 1);
  T * tmp = values;
  std::uninitialized_copy(&(vect[0]), &(vect[size-1]), tmp);

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
Array<T, is_scal>::~Array () {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG(dblAccessory, "Freeing "
	       << printMemorySize<T>(allocated_size*nb_component)
	       << " (" << id <<")");

  if(values){
    if(!is_scal)
      for (UInt i = 0; i < size * nb_component; ++i) {
	T * obj = values+i;
	obj->~T();
      }
    free(values);
  }
  size = allocated_size = 0;
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
void Array<T, is_scal>::allocate(UInt size,
                                 UInt nb_component) {
  AKANTU_DEBUG_IN();
  if (size == 0){
    values = NULL;
  } else {
    values = static_cast<T*>(malloc(nb_component * size * sizeof(T)));
    AKANTU_DEBUG_ASSERT(values != NULL,
			"Cannot allocate "
			<< printMemorySize<T>(size*nb_component)
			<< " (" << id <<")");
  }

  if (values == NULL) {
    this->size = this->allocated_size = 0;
  } else {
    AKANTU_DEBUG(dblAccessory, "Allocated "
		 << printMemorySize<T>(size*nb_component)
		 << " (" << id <<")");
    this->size = this->allocated_size = size;
  }

  this->size_of_type = sizeof(T);
  this->nb_component = nb_component;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * change the size of the array and allocate or free memory if needed. If the
 * size increases, the new tuples are filled with zeros
 * @param size new number of tuples contained in the array */
template <class T, bool is_scal>
void Array<T, is_scal>::resize(UInt new_size) {
  UInt old_size = size;

  T * old_values = values;
  if(new_size < size) {
    for (UInt i = new_size * nb_component; i < size * nb_component; ++i) {
      T * obj = old_values+i;
      obj->~T();
    }
  }

  resizeUnitialized(new_size);


  T val = T();
  if(size > old_size)
    std::uninitialized_fill(values + old_size*nb_component, values + size*nb_component, val);
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
void Array<T, is_scal>::resizeUnitialized(UInt new_size) {
  //  AKANTU_DEBUG_IN();
  // free some memory
  if(new_size <= allocated_size) {
    if(allocated_size - new_size > AKANTU_MIN_ALLOCATION) {
      AKANTU_DEBUG(dblAccessory, "Freeing "
		   << printMemorySize<T>((allocated_size - size)*nb_component)
		   << " (" << id <<")");

      // Normally there are no allocation problem when reducing an array
      T * tmp_ptr = static_cast<T*>(realloc(values, new_size * nb_component * sizeof(T)));
      if(new_size != 0 && tmp_ptr == NULL) {
	AKANTU_DEBUG_ERROR("Cannot free data (" << id << ")"
			   << " [current allocated size : " << allocated_size << " | "
			   << "requested size : " << new_size << "]");
      }
      values = tmp_ptr;
      allocated_size = new_size;
    }

    size = new_size;

    //    AKANTU_DEBUG_OUT();
    return;
  }

  // allocate more memory
  UInt size_to_alloc = (new_size - allocated_size < AKANTU_MIN_ALLOCATION) ?
    allocated_size + AKANTU_MIN_ALLOCATION : new_size;

  T *tmp_ptr = static_cast<T*>(realloc(values, size_to_alloc * nb_component * sizeof(T)));
  AKANTU_DEBUG_ASSERT(tmp_ptr != NULL,
                      "Cannot allocate "
		      << printMemorySize<T>(size_to_alloc * nb_component));
  if (tmp_ptr == NULL) {
    AKANTU_DEBUG_ERROR("Cannot allocate more data (" << id << ")"
		       << " [current allocated size : " << allocated_size << " | "
		       << "requested size : " << new_size << "]");
  }

  AKANTU_DEBUG(dblAccessory, "Allocating "
	       << printMemorySize<T>((size_to_alloc - allocated_size)*nb_component));

  allocated_size = size_to_alloc;
  size = new_size;
  values = tmp_ptr;

  //  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
void Array<T, is_scal>::extendComponentsInterlaced(UInt multiplicator,
                                                   UInt block_size) {
  AKANTU_DEBUG_IN();

  if (multiplicator == 1) return;

  AKANTU_DEBUG_ASSERT(multiplicator > 1,
  		      "invalid multiplicator");
  AKANTU_DEBUG_ASSERT(nb_component%block_size == 0,
		      "stride must divide actual number of components");

  values = static_cast<T*>(realloc(values, nb_component*multiplicator*size* sizeof(T)));

  UInt new_component = nb_component/block_size * multiplicator;

  for (UInt i = 0,k=size-1; i < size; ++i,--k) {
    for (UInt j = 0; j < new_component; ++j) {
      UInt m = new_component - j -1;
      UInt n = m/multiplicator;
      for (UInt l = 0,p=block_size-1;  l < block_size; ++l,--p) {
	values[k*nb_component*multiplicator+m*block_size+p] =
	  values[k*nb_component+n*block_size+p];
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
Int Array<T, is_scal>::find(const T & elem) const {
  AKANTU_DEBUG_IN();
  UInt i = 0;
  for (; (i < size) && (values[i] != elem); ++i);

  AKANTU_DEBUG_OUT();
  return (i == size) ? -1 : (Int) i;
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
Int Array<T, is_scal>::find(T elem[]) const {
  AKANTU_DEBUG_IN();
  T * it = values;
  UInt i = 0;
  for (;i < size; ++i) {
    if(*it == elem[0]) {
      T * cit = it;
      UInt c = 0;
      for(; (c < nb_component) && (*cit == elem[c]); ++c, ++cit);
      if(c == nb_component) {
	AKANTU_DEBUG_OUT();
	return i;
      }
    }
    it += nb_component;
  }
  return -1;
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
void Array<T, is_scal>::copy(const Array<T, is_scal>& vect, bool no_sanity_check) {
  AKANTU_DEBUG_IN();

  if(!no_sanity_check)
    if(vect.nb_component != nb_component)
      AKANTU_DEBUG_ERROR("The two arrays do not have the same number of components");

  resize((vect.size * vect.nb_component) / nb_component);

  T * tmp = values;
  std::uninitialized_copy(vect.storage(), vect.storage() + size * nb_component, tmp);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<bool is_scal>
class ArrayPrintHelper {
public:
  template<typename T>
  static void print_content(const Array<T> & vect, std::ostream & stream, int indent) {
    if(AKANTU_DEBUG_TEST(dblDump)) {
      std::string space;
      for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

      stream << space << " + values         : {";
      for (UInt i = 0; i < vect.getSize(); ++i) {
	stream << "{";
	for (UInt j = 0; j < vect.getNbComponent(); ++j) {
	  stream << vect(i, j);
	  if(j != vect.getNbComponent() - 1) stream << ", ";
	}
	stream << "}";
	if(i != vect.getSize() - 1) stream << ", ";
      }
      stream << "}" << std::endl;
    }
  }
};

template<>
class ArrayPrintHelper<false> {
public:
  template<typename T>
  static void print_content(__attribute__((unused)) const Array<T> & vect,
			    __attribute__((unused)) std::ostream & stream,
			    __attribute__((unused)) int indent) { }
};


/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
void Array<T, is_scal>::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  std::streamsize prec        = stream.precision();
  std::ios_base::fmtflags ff  = stream.flags();

  stream.setf (std::ios_base::showbase);
  stream.precision(2);

  stream << space << "Array<" << debug::demangle(typeid(T).name()) << "> [" << std::endl;
  stream << space << " + id             : " << this->id << std::endl;
  stream << space << " + size           : " << this->size << std::endl;
  stream << space << " + nb_component   : " << this->nb_component << std::endl;
  stream << space << " + allocated size : " << this->allocated_size << std::endl;
  stream << space << " + memory size    : "
	 << printMemorySize<T>(allocated_size*nb_component) << std::endl;
  if(!AKANTU_DEBUG_LEVEL_IS_TEST())
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

inline void ArrayBase::empty() {
  size = 0;
}

/* -------------------------------------------------------------------------- */
/* Iterators                                                                  */
/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
template<class R, class IR, bool is_r_scal>
class Array<T, is_scal>::iterator_internal {
public:
  typedef R                               value_type;
  typedef R*                              pointer;
  typedef R&                              reference;
  typedef const R&                        const_reference;
  typedef IR                              internal_value_type;
  typedef IR*                             internal_pointer;
  typedef std::ptrdiff_t                  difference_type;
  typedef std::random_access_iterator_tag iterator_category;

public:
  iterator_internal() : _offset(0), initial(NULL), ret(NULL), ret_ptr(NULL) {};

  iterator_internal(pointer_type data, UInt _offset)  :
    _offset(_offset),
    initial(data),
    ret(NULL),
    ret_ptr(data) {
    AKANTU_DEBUG_ERROR("The constructor should never be called it is just an ugly trick...");
  }

  iterator_internal(pointer wrapped)  : _offset(wrapped->size()),
					initial(wrapped->storage()),
					ret(const_cast<internal_pointer>(wrapped)),
					ret_ptr(wrapped->storage()) {
  }

  iterator_internal(const iterator_internal & it) {
    if(this != &it) {
      this->_offset = it._offset;
      this->initial = it.initial;
      this->ret_ptr = it.ret_ptr;
      this->ret = new internal_value_type(*it.ret, false);
    }
  }

  virtual ~iterator_internal() { delete ret; };

  inline iterator_internal & operator=(const iterator_internal & it) {
    if(this != &it) {
      this->_offset = it._offset;
      this->initial = it.initial;
      this->ret_ptr = it.ret_ptr;
      if(this->ret) this->ret->shallowCopy(*it.ret);
      else this->ret = new internal_value_type(*it.ret, false);
    }
    return *this;
  }

  UInt getCurrentIndex(){return (this->ret_ptr - this->initial)/this->_offset;};

  inline reference operator*() { ret->values = ret_ptr; return *ret; };
  inline const_reference operator*() const { ret->values = ret_ptr; return *ret; };
  inline pointer operator->() { ret->values = ret_ptr; return ret; };
  inline iterator_internal & operator++() { ret_ptr += _offset; return *this; };
  inline iterator_internal & operator--() { ret_ptr -= _offset; return *this; };

  inline iterator_internal & operator+=(const UInt n) { ret_ptr += _offset * n; return *this; }
  inline iterator_internal & operator-=(const UInt n) { ret_ptr -= _offset * n; return *this; }

  inline reference operator[](const UInt n) { ret->values = ret_ptr + n*_offset; return *ret; }
  inline const_reference operator[](const UInt n) const { ret->values = ret_ptr + n*_offset; return *ret; }

  inline bool operator==(const iterator_internal & other) const { return this->ret_ptr == other.ret_ptr; }
  inline bool operator!=(const iterator_internal & other) const { return this->ret_ptr != other.ret_ptr; }
  inline bool operator <(const iterator_internal & other) const { return this->ret_ptr  < other.ret_ptr; }
  inline bool operator<=(const iterator_internal & other) const { return this->ret_ptr <= other.ret_ptr; }
  inline bool operator> (const iterator_internal & other) const { return this->ret_ptr >  other.ret_ptr; }
  inline bool operator>=(const iterator_internal & other) const { return this->ret_ptr >= other.ret_ptr; }

  inline iterator_internal operator+(difference_type n) { iterator_internal tmp(*this); tmp += n; return tmp; }
  inline iterator_internal operator-(difference_type n) { iterator_internal tmp(*this); tmp -= n; return tmp; }


  inline difference_type operator-(const iterator_internal & b) { return (this->ret_ptr - b.ret_ptr) / _offset; }


  inline pointer_type data() const { return ret_ptr; }
  inline difference_type offset() const { return _offset; }

protected:
  UInt _offset;
  pointer_type initial;
  internal_pointer ret;
  pointer_type ret_ptr;
};

/* -------------------------------------------------------------------------- */
/**
 * Specialization for scalar types
 */
template <class T, bool is_scal>
template <class R, class IR>
class Array<T, is_scal>::iterator_internal<R, IR, true> {
public:
  typedef R                               value_type;
  typedef R*                              pointer;
  typedef R&                              reference;
  typedef const R&                        const_reference;
  typedef IR                              internal_value_type;
  typedef IR*                             internal_pointer;
  typedef std::ptrdiff_t                  difference_type;
  typedef std::random_access_iterator_tag iterator_category;

public:
  iterator_internal(pointer data = NULL, __attribute__ ((unused)) UInt _offset = 1) : _offset(_offset), ret(data), initial(data) { };
  iterator_internal(const iterator_internal & it) {
    if(this != &it) { this->ret = it.ret; this->initial = it.initial; }
  }

  virtual ~iterator_internal() { };

  inline iterator_internal & operator=(const iterator_internal & it)
  { if(this != &it) { this->ret = it.ret; this->initial = it.initial; } return *this; }

  UInt getCurrentIndex(){return (this->ret - this->initial)/this->_offset;};

  inline reference operator*() { return *ret; };
  inline const_reference operator*() const { return *ret; };
  inline pointer operator->() { return ret; };
  inline iterator_internal & operator++() { ++ret; return *this; };
  inline iterator_internal & operator--() { --ret; return *this; };

  inline iterator_internal & operator+=(const UInt n) { ret += n; return *this; }
  inline iterator_internal & operator-=(const UInt n) { ret -= n; return *this; }

  inline reference operator[](const UInt n) { return ret[n]; }

  inline bool operator==(const iterator_internal & other) const { return ret == other.ret; }
  inline bool operator!=(const iterator_internal & other) const { return ret != other.ret; }
  inline bool operator< (const iterator_internal & other) const { return ret <  other.ret; }
  inline bool operator<=(const iterator_internal & other) const { return ret <= other.ret; }
  inline bool operator> (const iterator_internal & other) const { return ret >  other.ret; }
  inline bool operator>=(const iterator_internal & other) const { return ret >= other.ret; }

  inline iterator_internal operator-(difference_type n) { return iterator_internal(ret - n); }
  inline iterator_internal operator+(difference_type n) { return iterator_internal(ret + n); }

  inline difference_type operator-(const iterator_internal & b) { return ret - b.ret; }

  inline pointer data() const { return ret; }
  inline difference_type offset() const { return _offset; }
protected:
  difference_type _offset;
  pointer ret;
  pointer initial;
};

/* -------------------------------------------------------------------------- */
/* Begin/End functions implementation                                         */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/**
 * Get an iterator that behaves like a pointer akantu::Vector<T> * to the
 * first tuple of the array.
 * @param n vector size. Has to be equal to nb_component. This unfortunate
 * redundancy is necessary to distinguish it from ::begin() which it
 * overloads. If compiled in debug mode, an incorrect value of n will result
 * in an exception being thrown. Optimized code will fail in an unpredicted
 * manner.
 * @return a vector_iterator
 */
template <class T, bool is_scal>
inline typename Array<T, is_scal>::vector_iterator Array<T, is_scal>::begin(UInt n) {
  AKANTU_DEBUG_ASSERT(nb_component == n,
		      "The iterator is not compatible with the type Array("
		      << n<< ")");
  return vector_iterator(new Vector<T>(values, n));
}

/* -------------------------------------------------------------------------- */
/**
 * Get an iterator that behaves like a pointer akantu::Vector<T> * pointing
 * *past* the last tuple of the array.
 * @param n vector size. see Array::begin(UInt n) for more
 * @return a vector_iterator
 */
template <class T, bool is_scal>
inline typename Array<T, is_scal>::vector_iterator Array<T, is_scal>::end(UInt n) {
  AKANTU_DEBUG_ASSERT(nb_component == n,
		      "The iterator is not compatible with the type Array("
		      << n<< ")");
  return vector_iterator(new Vector<T>(values + nb_component * size,
                                             n));
}

/* -------------------------------------------------------------------------- */
/**
 * Get a const iterator that behaves like a pointer akantu::Vector<T> * to the
 * first tuple of the array.
 * @param n vector size. see Array::begin(UInt n) for more
 * @return a vector_iterator
 */
template <class T, bool is_scal>
inline typename Array<T, is_scal>::const_vector_iterator Array<T, is_scal>::begin(UInt n) const {
  AKANTU_DEBUG_ASSERT(nb_component == n,
		      "The iterator is not compatible with the type Array("
		      << n<< ")");
  return const_vector_iterator(new Vector<T>(values, n));
}

/* -------------------------------------------------------------------------- */
/**
 * Get a const iterator that behaves like a pointer akantu::Vector<T> * pointing
 * *past* the last tuple of the array.
 * @param n vector size. see Array::begin(UInt n) for more
 * @return a const_vector_iterator
 */
template <class T, bool is_scal>
inline typename Array<T, is_scal>::const_vector_iterator Array<T, is_scal>::end(UInt n) const {
  AKANTU_DEBUG_ASSERT(nb_component == n,
		      "The iterator is not compatible with the type Array("
		      << n<< ")");
  return const_vector_iterator(new Vector<T>(values + nb_component * size,
                                                   n));
}

/* -------------------------------------------------------------------------- */
/**
 * Get an iterator that behaves like a pointer akantu::Vector<T> * to the
 * first tuple of the array.
 *
 * The reinterpret iterators allow to iterate over an array in any way that
 * preserves the number of entries of the array. This can for instance be use
 * full if the shape of the data in an array is not initially known.
 * @param n vector size.
 * @param size number of tuples in array. n times size must match the number
 * of entries of the array. If compiled in debug mode, an incorrect
 * combination of n and size will result
 * in an exception being thrown. Optimized code will fail in an unpredicted
 * manner.
 * @return a vector_iterator
 */
template <class T, bool is_scal>
inline typename Array<T, is_scal>::vector_iterator
Array<T, is_scal>::begin_reinterpret(UInt n, __attribute__((unused)) UInt size) {
  AKANTU_DEBUG_ASSERT(n * size == this->nb_component * this->size,
		      "The new values for size (" << size
		      << ") and nb_component (" << n
		      << ") are not compatible with the one of this array("
		      << this->size << "," << this->nb_component << ")");
  return vector_iterator(new Vector<T>(values, n));
}

/* -------------------------------------------------------------------------- */
/**
 * Get an iterator that behaves like a pointer akantu::Vector<T> * pointing
 * *past* the last tuple of the array.
 * @param n vector size.
 * @param size number of tuples in array. See Array::begin_reinterpret(UInt n, UInt size)
 * @return a vector_iterator
 */
template <class T, bool is_scal>
inline typename Array<T, is_scal>::vector_iterator
Array<T, is_scal>::end_reinterpret(UInt n, UInt size) {
  AKANTU_DEBUG_ASSERT(n * size == this->nb_component * this->size,
		      "The new values for size (" << size
		      << ") and nb_component (" << n
		      << ") are not compatible with the one of this array("
		      << this->size << "," << this->nb_component << ")");
  return vector_iterator(new Vector<T>(values + n * size, n));
}


/* -------------------------------------------------------------------------- */
/**
 * Get a const iterator that behaves like a pointer akantu::Vector<T> * to the
 * first tuple of the array.
 * @param n vector size.
 * @param size number of tuples in array. See Array::begin_reinterpret(UInt n, UInt size)
 * @return a const_vector_iterator
 */
template <class T, bool is_scal>
inline typename Array<T, is_scal>::const_vector_iterator
Array<T, is_scal>::begin_reinterpret(UInt n, __attribute__((unused)) UInt size) const {
  AKANTU_DEBUG_ASSERT(n * size == this->nb_component * this->size,
		      "The new values for size (" << size
		      << ") and nb_component (" << n
		      << ") are not compatible with the one of this array("
		      << this->size << "," << this->nb_component << ")");
  return const_vector_iterator(new Vector<T>(values, n));
}

/* -------------------------------------------------------------------------- */
/**
 * Get a const iterator that behaves like a pointer akantu::Vector<T> * pointing
 * *past* the last tuple of the array.
 * @param n vector size.
 * @param size number of tuples in array. See Array::begin_reinterpret(UInt n, UInt size)
 * @return a const_vector_iterator
 */
template <class T, bool is_scal>
inline typename Array<T, is_scal>::const_vector_iterator
Array<T, is_scal>::end_reinterpret(UInt n, UInt size) const {
  AKANTU_DEBUG_ASSERT(n * size == this->nb_component * this->size,
		      "The new values for size (" << size
		      << ") and nb_component (" << n
		      << ") are not compatible with the one of this array("
		      << this->size << "," << this->nb_component << ")");
  return const_vector_iterator(new Vector<T>(values + n * size, n));
}

/* -------------------------------------------------------------------------- */
/**
 * Get an iterator that behaves like a pointer akantu::Matrix<T> * to the
 * first tuple of the array.
 * @param m number of rows
 * @param n number of columns. m times n has to equal nb_component.
 * If compiled in debug mode, an incorrect combination of m and n will result
 * in an exception being thrown. Optimized code will fail in an unpredicted
 * manner.
 * @return a matrix_iterator
 */
template <class T, bool is_scal>
inline typename Array<T, is_scal>::matrix_iterator Array<T, is_scal>::begin(UInt m, UInt n) {
  AKANTU_DEBUG_ASSERT(nb_component == n*m,
		      "The iterator is not compatible with the type Matrix("
		      << m << "," << n<< ")");
  return matrix_iterator(new Matrix<T>(values, m, n));
}

/* -------------------------------------------------------------------------- */
/**
 * Get an iterator that behaves like a pointer akantu::Matrix<T> * pointing
 * *past* the last tuple of the array.
 * @param m number of rows
 * @param n number of columns. See Array::begin(UInt m, UInt n)
 * @return a matrix_iterator
 */
template <class T, bool is_scal>
inline typename Array<T, is_scal>::matrix_iterator Array<T, is_scal>::end(UInt m, UInt n) {
  AKANTU_DEBUG_ASSERT(nb_component == n*m,
		      "The iterator is not compatible with the type Matrix("
		      << m << "," << n<< ")");
  return matrix_iterator(new Matrix<T>(values + nb_component * size, m, n));
}

/* -------------------------------------------------------------------------- */
/**
 * Get a const iterator that behaves like a pointer akantu::Matrix<T> * to the
 * first tuple of the array.
 * @param m number of rows
 * @param n number of columns. See Array::begin(UInt m, UInt n)
 * @return a matrix_iterator
 */
template <class T, bool is_scal>
inline typename Array<T, is_scal>::const_matrix_iterator Array<T, is_scal>::begin(UInt m, UInt n) const {
  AKANTU_DEBUG_ASSERT(nb_component == n*m,
		      "The iterator is not compatible with the type Matrix("
		      << m << "," << n<< ")");
  return const_matrix_iterator(new Matrix<T>(values, m, n));
}

/* -------------------------------------------------------------------------- */
/**
 * Get a const iterator that behaves like a pointer akantu::Matrix<T> * pointing
 * *past* the last tuple of the array.
 * @param m number of rows
 * @param n number of columns. See Array::begin(UInt m, UInt n)
 * @return a const_matrix_iterator
 */
template <class T, bool is_scal>
inline typename Array<T, is_scal>::const_matrix_iterator Array<T, is_scal>::end(UInt m, UInt n) const {
  AKANTU_DEBUG_ASSERT(nb_component == n*m,
		      "The iterator is not compatible with the type Matrix("
		      << m << "," << n<< ")");
  return const_matrix_iterator(new Matrix<T>(values + nb_component * size, m, n));
}

/* -------------------------------------------------------------------------- */
/**
 * Get an iterator that behaves like a pointer akantu::Matrix<T> * to the
 * first tuple of the array.
 *
 * The reinterpret iterators allow to iterate over an array in any way that
 * preserves the number of entries of the array. This can for instance be use
 * full if the shape of the data in an array is not initially known.
 * @param m number of rows
 * @param n number of columns
 * @param size number of tuples in array. m times n times size must match the number
 * of entries of the array. If compiled in debug mode, an incorrect
 * combination of m, n and size will result
 * in an exception being thrown. Optimized code will fail in an unpredicted
 * manner.
 * @return a matrix_iterator
 */
template <class T, bool is_scal>
inline typename Array<T, is_scal>::matrix_iterator
Array<T, is_scal>::begin_reinterpret(UInt m, UInt n, __attribute__((unused)) UInt size) {
  AKANTU_DEBUG_ASSERT(n * m * size == this->nb_component * this->size,
		      "The new values for size (" << size
		      << ") and nb_component (" << m << "," << n << " = " << n * m
		      << ") are not compatible with the one of this array("
		      << this->size << "," << this->nb_component << ")");
  return matrix_iterator(new Matrix<T>(values, m, n));
}

/* -------------------------------------------------------------------------- */
/**
 * Get an iterator that behaves like a pointer akantu::Matrix<T> * pointing
 * *past* the last tuple of the array.
 * @param m number of rows
 * @param n number of columns
 * @param size number of tuples in array. See Array::begin_reinterpret(UInt m, UInt n, UInt size)
 * @return a matrix_iterator
 */
template <class T, bool is_scal>
inline typename Array<T, is_scal>::matrix_iterator
Array<T, is_scal>::end_reinterpret(UInt m, UInt n, UInt size) {
  AKANTU_DEBUG_ASSERT(n * m * size == this->nb_component * this->size,
		      "The new values for size (" << size
		      << ") and nb_component (" << m << "," << n << " = " << n * m
		      << ") are not compatible with the one of this array("
		      << this->size << "," << this->nb_component << ")");
  return matrix_iterator(new Matrix<T>(values + n * m * size, m, n));
}


/* -------------------------------------------------------------------------- */
/**
 * Get a const iterator that behaves like a pointer akantu::Matrix<T> * to the
 * first tuple of the array.
 * @param m number of rows
 * @param n number of columns
 * @param size number of tuples in array. See Array::begin_reinterpret(UInt m, UInt n, UInt size)
 * @return a const_matrix_iterator
 */
template <class T, bool is_scal>
inline typename Array<T, is_scal>::const_matrix_iterator
Array<T, is_scal>::begin_reinterpret(UInt m, UInt n, __attribute__((unused)) UInt size) const {
  AKANTU_DEBUG_ASSERT(n * m * size == this->nb_component * this->size,
		      "The new values for size (" << size
		      << ") and nb_component (" << m << "," << n << " = " << n * m
		      << ") are not compatible with the one of this array("
		      << this->size << "," << this->nb_component << ")");
  return const_matrix_iterator(new Matrix<T>(values, m, n));
}

/* -------------------------------------------------------------------------- */
/**
 * Get a const iterator that behaves like a pointer akantu::Matrix<T> * pointing
 * *past* the last tuple of the array.
 * @param m number of rows
 * @param n number of columns
 * @param size number of tuples in array. See Array::begin_reinterpret(UInt m, UInt n, UInt size)
 * @return a const_matrix_iterator
 */
template <class T, bool is_scal>
inline typename Array<T, is_scal>::const_matrix_iterator
Array<T, is_scal>::end_reinterpret(UInt m, UInt n, UInt size) const {
  AKANTU_DEBUG_ASSERT(n * m * size == this->nb_component * this->size,
		      "The new values for size (" << size
		      << ") and nb_component (" << m << "," << n << " = " << n * m
		      << ") are not compatible with the one of this array("
		      << this->size << "," << this->nb_component << ")");
  return const_matrix_iterator(new Matrix<T>(values + n * m * size, m, n));
}

/* -------------------------------------------------------------------------- */
/** Get an iterator that behaves like a pointer T * to the
 *  first entry in the member array values
 *  @return a scalar_iterator
 */
template <class T, bool is_scal>
inline Array<T, is_scal>::iterator<T> Array<T, is_scal>::begin() {
  AKANTU_DEBUG_ASSERT(nb_component == 1, "this iterator cannot be used on a vector which has nb_component != 1");
  return iterator<T>(values);
}

/* -------------------------------------------------------------------------- */
  /*! Get an iterator that behaves like a pointer T * that points *past* the
   *  last entry in the member array values
   *  @return a scalar_iterator
   */
template <class T, bool is_scal>
inline Array<T, is_scal>::iterator<T> Array<T, is_scal>::end() {
  AKANTU_DEBUG_ASSERT(nb_component == 1, "this iterator cannot be used on a array which has nb_component != 1");
  return iterator<T>(values + size);
}

/* -------------------------------------------------------------------------- */
/*! Get a const iterator that behaves like a pointer T * to the
 *  first entry in the member array values
 *  @return a const_scalar_iterator
 */
template <class T, bool is_scal>
inline Array<T, is_scal>::const_iterator<T> Array<T, is_scal>::begin() const {
  AKANTU_DEBUG_ASSERT(nb_component == 1, "this iterator cannot be used on a array which has nb_component != 1");
  return const_iterator<T>(values);
}

/* -------------------------------------------------------------------------- */
/*! Get a const iterator that behaves like a pointer T * that points *past* the
 *  last entry in the member array values
 *  @return a const_scalar_iterator
 */
template <class T, bool is_scal>
inline Array<T, is_scal>::const_iterator<T> Array<T, is_scal>::end() const {
  AKANTU_DEBUG_ASSERT(nb_component == 1, "this iterator cannot be used on a array which has nb_component != 1");
  return const_iterator<T>(values + size);
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
template<typename R>
class Array<T, is_scal>::const_iterator : public iterator_internal<const R, R> {
public:
  typedef iterator_internal<const R, R> parent;
  typedef typename parent::value_type	       value_type;
  typedef typename parent::pointer	       pointer;
  typedef typename parent::reference	       reference;
  typedef typename parent::difference_type   difference_type;
  typedef typename parent::iterator_category iterator_category;
public:

  const_iterator() : parent() {};
  const_iterator(pointer_type data, UInt offset) : parent(data, offset) {}
  const_iterator(pointer warped) : parent(warped) {}
  const_iterator(const parent & it) : parent(it) {}
  //    const_iterator(const const_iterator<R> & it) : parent(it) {}

  inline const_iterator operator+(difference_type n)
  { return parent::operator+(n); }
  inline const_iterator operator-(difference_type n)
  { return parent::operator-(n); }
  inline difference_type operator-(const const_iterator & b)
  { return parent::operator-(b); }

  inline const_iterator & operator++()
  { parent::operator++(); return *this; };
  inline const_iterator & operator--()
  { parent::operator--(); return *this; };
  inline const_iterator & operator+=(const UInt n)
  { parent::operator+=(n); return *this; }
};
// #endif


// #if defined(AKANTU_CORE_CXX11)
//   template<class R> using iterator = iterator_internal<R>;
// #else
template < class T, class R, bool issame = is_same<T, R>::value >
struct ConstConverterIteratorHelper {
  typedef typename Array<T>::template const_iterator<R> const_iterator;
  typedef typename Array<T>::template iterator<R> iterator;
  static inline const_iterator convert(const iterator & it) {
    return const_iterator(new R(*it, false));
  }
};

template < class T, class R >
struct ConstConverterIteratorHelper<T, R, true> {
  typedef typename Array<T>::template const_iterator<R> const_iterator;
  typedef typename Array<T>::template iterator<R> iterator;
  static inline const_iterator convert(const iterator & it) {
    return const_iterator(it.data(), it.offset());
  }
};


template <class T, bool is_scal>
template<typename R>
class  Array<T, is_scal>::iterator : public iterator_internal<R> {
public:
  typedef iterator_internal<R> parent;
  typedef typename parent::value_type	       value_type;
  typedef typename parent::pointer	       pointer;
  typedef typename parent::reference	       reference;
  typedef typename parent::difference_type   difference_type;
  typedef typename parent::iterator_category iterator_category;
public:
  iterator() : parent() {};
  iterator(pointer_type data, UInt offset) : parent(data, offset) {};
  iterator(pointer warped) : parent(warped) {}
  iterator(const parent & it) : parent(it) {}
  //    iterator(const iterator<R> & it) : parent(it) {}

  operator const_iterator<R>() {
    return ConstConverterIteratorHelper<T, R>::convert(*this);
  }

  inline iterator operator+(difference_type n)
  { return parent::operator+(n);; }
  inline iterator operator-(difference_type n)
  { return parent::operator-(n);; }
  inline difference_type operator-(const iterator & b)
  { return parent::operator-(b); }

  inline iterator & operator++()
  { parent::operator++(); return *this; };
  inline iterator & operator--()
  { parent::operator--(); return *this; };
  inline iterator & operator+=(const UInt n)
  { parent::operator+=(n); return *this; }
};

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
template<typename R>
inline Array<T, is_scal>::iterator<R> Array<T, is_scal>::erase(const iterator<R> & it) {
  T * curr = it.data();
  UInt pos = (curr - values) / nb_component;
  erase(pos);
  iterator<R> rit = it;
  return --rit;
}


// #endif

