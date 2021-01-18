/**
 * @file   synchronized_array.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Tue Dec 02 2014
 * @date last modification: Fri Feb 23 2018
 *
 * @brief  synchronized array: a array can be registered to another (hereafter
 * called top) array. If an element is added to or removed from the top array,
 * the registered array removes or adds at the same position an element. The two
 * arrays stay therefore synchronized.
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#ifndef AST_SYNCHRONIZED_ARRAY_HH_
#define AST_SYNCHRONIZED_ARRAY_HH_

/* -------------------------------------------------------------------------- */
// std
#include <vector>

// akantu
#include "aka_array.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
enum SyncChoice { _added, _deleted };

/* -------------------------------------------------------------------------- */
class SynchronizedArrayBase {
public:
  SynchronizedArrayBase() = default;
  ~SynchronizedArrayBase() = default;

  virtual ID getID() const { return "call should be virtual"; };

  virtual UInt syncDeletedElements(std::vector<UInt> & delete_el) = 0;
  virtual UInt syncAddedElements(UInt nb_added_el) = 0;
};

/* -------------------------------------------------------------------------- */
template <class T>
class SynchronizedArray : public SynchronizedArrayBase, protected Array<T> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  using value_type = typename Array<T>::value_type;
  using reference = typename Array<T>::reference;
  using pointer_type = typename Array<T>::pointer_type;
  using const_reference = typename Array<T>::const_reference;

  SynchronizedArray(UInt size, UInt nb_component, const_reference value,
                    const ID & id, const_reference default_value,
                    const std::string & restart_name);
  ~SynchronizedArray() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// push_back
  template <typename P> inline void push_back(P && value);

  /// erase
  inline void erase(UInt i);
  // template<typename R>
  // inline void erase(const iterator<R> & it);

  /// synchronize elements
  void syncElements(SyncChoice sync_choice);

  /// dump restart file
  void dumpRestartFile(const std::string & file_name) const;

  /// read restart file
  void readRestartFile(const std::string & file_name);

  /// register depending array
  void registerDependingArray(SynchronizedArrayBase & array);

  /// function to print the contain of the class
  void printself(std::ostream & stream, int indent = 0) const override;

  /// find position of element
  Int find(const T & elem) const { return Array<T>::find(elem); };

  /// set values to zero
  inline void zero() { Array<T>::zero(); };
  //  inline void clear() { memset(values, 0, size*nb_component*sizeof(T)); };

  /// set all entries of the array to the value t
  /// @param t value to fill the array with
  inline void set(T t) { Array<T>::set(t); }

  /// set
  template <template <typename> class C> inline void set(const C<T> & vm) {
    Array<T>::set(vm);
  }

  /// set all entries of the array to value t and set default value
  inline void setAndChangeDefault(T t) {
    this->set(t);
    this->default_value = t;
  }

  /// copy the content of an other array
  void copy(const SynchronizedArray<T> & vect) { Array<T>::copy(vect); };

  /// give the address of the memory allocated for this array
  T * storage() const { return Array<T>::storage(); };
  //  T * storage() const { return this->values; };

  // get nb component
  UInt getNbComponent() const { return Array<T>::getNbComponent(); };

protected:
  UInt syncDeletedElements(std::vector<UInt> & del_elements) override;
  UInt syncAddedElements(UInt nb_add_elements) override;

  /* ------------------------------------------------------------------------ */
  /* Operators                                                                */
  /* ------------------------------------------------------------------------ */
public:
  inline reference operator()(UInt i, UInt j = 0) {
    return Array<T>::operator()(i, j);
  }
  inline const_reference operator()(UInt i, UInt j = 0) const {
    return Array<T>::operator()(i, j);
  }

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_SET_MACRO(DefaultValue, default_value, T);

  UInt size() const { return this->size_; };

  ID getID() const override { return Array<T>::getID(); };

  const Array<T> & getArray() const {
    const Array<T> & a = *(dynamic_cast<const Array<T> *>(this));
    return a;
  };

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// init value when new elements added
  T default_value;

  /// restart file_name
  const std::string restart_name;

  /// elements that have been deleted
  std::vector<UInt> deleted_elements;

  /// number of elements to add
  UInt nb_added_elements;

  /// pointers to arrays to be updated
  std::vector<SynchronizedArrayBase *> depending_arrays;
};

/// standard output stream operator
template <typename T>
inline std::ostream & operator<<(std::ostream & stream,
                                 const SynchronizedArray<T> & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "synchronized_array_inline_impl.hh"

#endif /* AST_SYNCHRONIZED_ARRAY_HH_ */
