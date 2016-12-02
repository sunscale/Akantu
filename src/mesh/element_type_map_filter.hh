/**
 * @file   element_type_map_filter.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Sep 02 2014
 * @date last modification: Fri Dec 18 2015
 *
 * @brief  Filtered version based on a an akantu::ElementGroup of a
 * akantu::ElementTypeMap
 *
 * @section LICENSE
 *
 * Copyright  (©)  2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de Lausanne)
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

#ifndef __AKANTU_BY_ELEMENT_TYPE_FILTER_HH__
#define __AKANTU_BY_ELEMENT_TYPE_FILTER_HH__
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/* ArrayFilter                                                                */
/* -------------------------------------------------------------------------- */

template <typename T> class ArrayFilter {

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */

public:
  /// standard iterator
  template <typename R = T> class iterator {
    inline bool operator!=(__attribute__((unused)) iterator<R> & other) { throw; };
    inline bool operator==(__attribute__((unused)) iterator<R> & other) { throw; };

    inline iterator<R> & operator++() { throw; };
    inline T operator*() {
      throw;
      return T();
    };
  };

  /// const iterator
  template <template <class S> class original_iterator, typename Shape,
            typename filter_iterator>
  class const_iterator {

  public:
    UInt getCurrentIndex() {
      return (*this->filter_it * this->nb_item_per_elem +
              this->sub_element_counter);
    }

    inline const_iterator(){};
    inline const_iterator(const original_iterator<Shape> & origin_it,
                          const filter_iterator & filter_it,
                          UInt nb_item_per_elem)
        : origin_it(origin_it), filter_it(filter_it),
          nb_item_per_elem(nb_item_per_elem), sub_element_counter(0){};

    inline bool operator!=(const_iterator & other) const {
      return !((*this) == other);
    }
    inline bool operator==(const_iterator & other) const {
      return (this->origin_it == other.origin_it &&
              this->filter_it == other.filter_it &&
              this->sub_element_counter == other.sub_element_counter);
    }

    inline bool operator!=(const const_iterator & other) const {
      return !((*this) == other);
    }
    inline bool operator==(const const_iterator & other) const {
      return (this->origin_it == other.origin_it &&
              this->filter_it == other.filter_it &&
              this->sub_element_counter == other.sub_element_counter);
    }

    inline const_iterator & operator++() {

      ++sub_element_counter;
      if (sub_element_counter == nb_item_per_elem) {
        sub_element_counter = 0;
        ++filter_it;
      }
      return *this;
    };

    inline Shape operator*() {
      return origin_it[nb_item_per_elem * (*filter_it) + sub_element_counter];
    };

  private:
    original_iterator<Shape> origin_it;
    filter_iterator filter_it;

    /// the number of item per element
    UInt nb_item_per_elem;
    /// counter for every sub element group
    UInt sub_element_counter;
  };

  typedef iterator<Vector<T> > vector_iterator;

  typedef Array<T> array_type;

  typedef const_iterator<array_type::template const_iterator, Vector<T>,
                         Array<UInt>::const_iterator<UInt> >
      const_vector_iterator;

  typedef typename array_type::value_type value_type;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

public:
  ArrayFilter(const Array<T> & array, const Array<UInt> & filter,
              UInt nb_item_per_elem)
      : array(array), filter(filter), nb_item_per_elem(nb_item_per_elem){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
  const_vector_iterator begin_reinterpret(UInt n, UInt new_size) const {
    AKANTU_DEBUG_ASSERT(
        n * new_size == this->getNbComponent() * this->getSize(),
        "The new values for size ("
            << new_size << ") and nb_component (" << n
            << ") are not compatible with the one of this array("
            << this->getSize() << "," << this->getNbComponent() << ")");
    UInt new_full_array_size =
        this->array.getSize() * array.getNbComponent() / n;
    UInt new_nb_item_per_elem = this->nb_item_per_elem;
    if (new_size != 0 && n != 0)
      new_nb_item_per_elem = this->array.getNbComponent() *
                             this->filter.getSize() * this->nb_item_per_elem /
                             (n * new_size);

    return const_vector_iterator(
        this->array.begin_reinterpret(n, new_full_array_size),
        this->filter.begin(), new_nb_item_per_elem);
  };

  const_vector_iterator end_reinterpret(UInt n, UInt new_size) const {
    AKANTU_DEBUG_ASSERT(
        n * new_size == this->getNbComponent() * this->getSize(),
        "The new values for size ("
            << new_size << ") and nb_component (" << n
            << ") are not compatible with the one of this array("
            << this->getSize() << "," << this->getNbComponent() << ")");
    UInt new_full_array_size =
        this->array.getSize() * this->array.getNbComponent() / n;
    UInt new_nb_item_per_elem = this->nb_item_per_elem;
    if (new_size != 0 && n != 0)
      new_nb_item_per_elem = this->array.getNbComponent() *
                             this->filter.getSize() * this->nb_item_per_elem /
                             (n * new_size);

    return const_vector_iterator(
        this->array.begin_reinterpret(n, new_full_array_size),
        this->filter.end(), new_nb_item_per_elem);
  };

  vector_iterator begin_reinterpret(UInt, UInt) { throw; };

  vector_iterator end_reinterpret(UInt, UInt) { throw; };

  /// return the size of the filtered array which is the filter size
  UInt getSize() const {
    return this->filter.getSize() * this->nb_item_per_elem;
  };
  /// the number of components of the filtered array
  UInt getNbComponent() const { return this->array.getNbComponent(); };

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

private:
  /// reference to array of data
  const Array<T> & array;
  /// reference to the filter used to select elements
  const Array<UInt> & filter;
  /// the number of item per element
  UInt nb_item_per_elem;
};

/* -------------------------------------------------------------------------- */
/* ElementTypeMapFilter */
/* -------------------------------------------------------------------------- */

template <class T, typename SupportType = ElementType>
class ElementTypeMapArrayFilter {

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */

public:
  typedef T type;
  typedef ArrayFilter<T> array_type;
  typedef typename array_type::value_type value_type;

  typedef typename ElementTypeMapArray<UInt, SupportType>::type_iterator
      type_iterator;

  // class type_iterator{

  // public:

  //    typedef typename ElementTypeMapArray<T,SupportType>::type_iterator
  //    type_it;

  //  public:
  //   type_iterator(){};
  //   //    type_iterator(const type_iterator & it){original_it =
  //   it.original_it;};
  //   type_iterator(const type_it & it){original_it = it;};

  //   inline ElementType & operator*(){throw;};
  //   inline ElementType & operator*() const{throw;};
  //   inline type_iterator & operator++(){throw;return *this;};
  //   type_iterator operator++(int){throw; return *this;};
  //   inline bool operator==(const type_iterator & other) const{throw;return
  //   false;};
  //   inline bool operator!=(const type_iterator & other) const{throw;return
  //   false;};
  //   //    type_iterator & operator=(const type_iterator & other){throw;return
  //   *this;};

  //   type_it original_it;

  // };

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

public:
  ElementTypeMapArrayFilter(
      const ElementTypeMapArray<T, SupportType> & array,
      const ElementTypeMapArray<UInt, SupportType> & filter,
      const ElementTypeMap<UInt, SupportType> & nb_data_per_elem)
      : array(array), filter(filter), nb_data_per_elem(nb_data_per_elem) {}

  ElementTypeMapArrayFilter(
      const ElementTypeMapArray<T, SupportType> & array,
      const ElementTypeMapArray<UInt, SupportType> & filter)
      : array(array), filter(filter) {}

  ~ElementTypeMapArrayFilter() {}

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

  inline const ArrayFilter<T>
  operator()(const SupportType & type,
             const GhostType & ghost_type = _not_ghost) const {
    if (filter.exists(type, ghost_type)) {
      if (nb_data_per_elem.exists(type, ghost_type))
        return ArrayFilter<T>(array(type, ghost_type), filter(type, ghost_type),
                              nb_data_per_elem(type, ghost_type) /
                                  array(type, ghost_type).getNbComponent());
      else
        return ArrayFilter<T>(array(type, ghost_type), filter(type, ghost_type),
                              1);
    } else {
      return ArrayFilter<T>(empty_array, empty_filter, 1);
    }
  };

  inline type_iterator firstType(UInt dim = _all_dimensions,
                                 GhostType ghost_type = _not_ghost,
                                 ElementKind kind = _ek_not_defined) const {
    return filter.firstType(dim, ghost_type, kind);
  };

  inline type_iterator lastType(UInt dim = _all_dimensions,
                                GhostType ghost_type = _not_ghost,
                                ElementKind kind = _ek_not_defined) const {
    return filter.lastType(dim, ghost_type, kind);
  };

  ElementTypeMap<UInt>
  getNbComponents(UInt dim = _all_dimensions, GhostType ghost_type = _not_ghost,
                  ElementKind kind = _ek_not_defined) const {
    return this->array.getNbComponents(dim, ghost_type, kind);
  };

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */

  std::string getID() {
    return std::string("filtered:" + this->array().getID());
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

protected:
  const ElementTypeMapArray<T, SupportType> & array;
  const ElementTypeMapArray<UInt, SupportType> & filter;
  ElementTypeMap<UInt> nb_data_per_elem;

  /// Empty array to be able to return consistent filtered arrays
  Array<T> empty_array;
  Array<UInt> empty_filter;
};

__END_AKANTU__

#endif /* __AKANTU_BY_ELEMENT_TYPE_FILTER_HH__ */
