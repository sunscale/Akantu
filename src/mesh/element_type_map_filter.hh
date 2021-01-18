/**
 * @file   element_type_map_filter.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Sep 02 2014
 * @date last modification: Sun Dec 03 2017
 *
 * @brief  Filtered version based on a an akantu::ElementGroup of a
 * akantu::ElementTypeMap
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef AKANTU_BY_ELEMENT_TYPE_FILTER_HH_
#define AKANTU_BY_ELEMENT_TYPE_FILTER_HH_
/* -------------------------------------------------------------------------- */

namespace akantu {

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
    inline bool operator!=(__attribute__((unused)) iterator<R> & other) {
      throw;
    };
    inline bool operator==(__attribute__((unused)) iterator<R> & other) {
      throw;
    };

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

    inline const_iterator() = default;
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

  using vector_iterator = iterator<Vector<T>>;

  using array_type = Array<T>;

  using const_vector_iterator =
      const_iterator<array_type::template const_iterator, Vector<T>,
                     Array<UInt>::const_iterator<UInt>>;

  using value_type = typename array_type::value_type;

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
        n * new_size == this->getNbComponent() * this->size(),
        "The new values for size ("
            << new_size << ") and nb_component (" << n
            << ") are not compatible with the one of this array("
            << this->size() << "," << this->getNbComponent() << ")");
    UInt new_full_array_size = this->array.size() * array.getNbComponent() / n;
    UInt new_nb_item_per_elem = this->nb_item_per_elem;
    if (new_size != 0 && n != 0) {
      new_nb_item_per_elem = this->array.getNbComponent() *
                             this->filter.size() * this->nb_item_per_elem /
                             (n * new_size);
    }

    return const_vector_iterator(
        this->array.begin_reinterpret(n, new_full_array_size),
        this->filter.begin(), new_nb_item_per_elem);
  };

  const_vector_iterator end_reinterpret(UInt n, UInt new_size) const {
    AKANTU_DEBUG_ASSERT(
        n * new_size == this->getNbComponent() * this->size(),
        "The new values for size ("
            << new_size << ") and nb_component (" << n
            << ") are not compatible with the one of this array("
            << this->size() << "," << this->getNbComponent() << ")");
    UInt new_full_array_size =
        this->array.size() * this->array.getNbComponent() / n;
    UInt new_nb_item_per_elem = this->nb_item_per_elem;
    if (new_size != 0 && n != 0) {
      new_nb_item_per_elem = this->array.getNbComponent() *
                             this->filter.size() * this->nb_item_per_elem /
                             (n * new_size);
    }

    return const_vector_iterator(
        this->array.begin_reinterpret(n, new_full_array_size),
        this->filter.end(), new_nb_item_per_elem);
  };

  vector_iterator begin_reinterpret(UInt /*unused*/, UInt /*unused*/) {
    throw;
  };

  vector_iterator end_reinterpret(UInt /*unused*/, UInt /*unused*/) { throw; };

  /// return the size of the filtered array which is the filter size
  UInt size() const { return this->filter.size() * this->nb_item_per_elem; };
  /// the number of components of the filtered array
  UInt getNbComponent() const { return this->array.getNbComponent(); };

  bool empty() const __attribute__((warn_unused_result)) { return (size() == 0);}

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
  using type = T;
  using array_type = ArrayFilter<T>;
  using value_type = typename array_type::value_type;

  using type_iterator =
      typename ElementTypeMapArray<UInt, SupportType>::type_iterator;

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

  ~ElementTypeMapArrayFilter() = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

  inline ArrayFilter<T> operator()(const SupportType & type,
                                   GhostType ghost_type = _not_ghost) const {
    if (filter.exists(type, ghost_type)) {
      if (nb_data_per_elem.exists(type, ghost_type)) {
        return ArrayFilter<T>(array(type, ghost_type), filter(type, ghost_type),
                              nb_data_per_elem(type, ghost_type) /
                                  array(type, ghost_type).getNbComponent());
      }
      return ArrayFilter<T>(array(type, ghost_type), filter(type, ghost_type),
                            1);
    }
    return ArrayFilter<T>(empty_array, empty_filter, 1);
  };

  template <typename... Args>
  decltype(auto) elementTypes(Args &&... args) const {
    return filter.elementTypes(std::forward<decltype(args)>(args)...);
  }

  decltype(auto) getNbComponents(UInt dim = _all_dimensions,
                                 GhostType ghost_type = _not_ghost,
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

} // namespace akantu

#endif /* AKANTU_BY_ELEMENT_TYPE_FILTER_HH_ */
