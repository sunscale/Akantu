/**
 * @file   dumper_element_iterator.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Sep 02 2014
 * @date last modification: Fri May 15 2015
 *
 * @brief  Iterators for elemental fields
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

#ifndef __AKANTU_DUMPER_ELEMENT_ITERATOR_HH__
#define __AKANTU_DUMPER_ELEMENT_ITERATOR_HH__
/* -------------------------------------------------------------------------- */
#include "element.hh"
/* -------------------------------------------------------------------------- */
namespace akantu {
__BEGIN_AKANTU_DUMPER__
/* -------------------------------------------------------------------------- */

template<class types, template <class> class final_iterator>
class element_iterator {
  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */
public:

  typedef typename types::it_type      it_type;
  typedef typename types::field_type field_type;
  typedef typename types::array_type array_type;
  typedef typename types::array_iterator array_iterator;
  typedef final_iterator<types> iterator;

public:

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

  element_iterator(const field_type & field,
                   const typename field_type::type_iterator & t_it,
                   const typename field_type::type_iterator & t_it_end,
                   const array_iterator & array_it,
                   const array_iterator & array_it_end,
                   const GhostType ghost_type = _not_ghost)
    : field(field),
      tit(t_it),
      tit_end(t_it_end),
      array_it(array_it),
      array_it_end(array_it_end),
      ghost_type(ghost_type) {
  }

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

public:

  bool operator!=(const iterator & it) const {
    return (ghost_type != it.ghost_type)
      || (tit != it.tit || (array_it != it.array_it));
  }

  iterator & operator++() {
    ++array_it;
    while(array_it == array_it_end && tit != tit_end) {
      ++tit;
      if(tit != tit_end) {

        const array_type & vect = field(*tit, ghost_type);
        UInt _nb_data_per_elem = getNbDataPerElem(*tit);
        UInt nb_component = vect.getNbComponent();
        UInt size = (vect.size() * nb_component) / _nb_data_per_elem;

        array_it       = vect.begin_reinterpret(_nb_data_per_elem,size);
        array_it_end   = vect.end_reinterpret  (_nb_data_per_elem,size);
      }
    }
    return *(static_cast<iterator *>(this));
  };

  ElementType getType() { return *tit; }

  UInt element_type() { return getIOHelperType(*tit); }

  Element getCurrentElement(){
    return Element(*tit,array_it.getCurrentIndex());
  }

  UInt getNbDataPerElem(const ElementType & type) const {
    if (!nb_data_per_elem.exists(type, ghost_type))
      return field(type,ghost_type).getNbComponent();
    
    return nb_data_per_elem(type,ghost_type);
  }

  void setNbDataPerElem(const ElementTypeMap<UInt> & nb_data){
    this->nb_data_per_elem = nb_data;
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

protected:

  /// the field to iterate on
  const field_type & field;
  /// field iterator
  typename field_type::type_iterator tit;
  /// field iterator end
  typename field_type::type_iterator tit_end;
  /// array iterator
  array_iterator array_it;
  /// internal iterator end
  array_iterator array_it_end;
  /// ghost type identification
  const GhostType ghost_type;
  /// number of data per element
  ElementTypeMap<UInt> nb_data_per_elem;

};

/* -------------------------------------------------------------------------- */
template<typename types>
class elemental_field_iterator
  : public element_iterator<types, elemental_field_iterator> {
public:
  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */

  typedef element_iterator<types, ::akantu::dumper::elemental_field_iterator> parent;
  typedef typename types::it_type     it_type;
  typedef typename types::return_type return_type;
  typedef typename types::field_type  field_type;
  typedef typename types::array_iterator array_iterator;

public:

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

  elemental_field_iterator(const field_type & field,
                           const typename field_type::type_iterator & t_it,
                           const typename field_type::type_iterator & t_it_end,
                           const array_iterator & array_it,
                           const array_iterator & array_it_end,
                           const GhostType ghost_type = _not_ghost) :
    parent(field, t_it, t_it_end, array_it, array_it_end, ghost_type) { }

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

  return_type operator*(){
    return *this->array_it;
  }

private:

};

/* -------------------------------------------------------------------------- */
__END_AKANTU_DUMPER__
} // akantu
/* -------------------------------------------------------------------------- */



#endif /* __AKANTU_DUMPER_ELEMENT_ITERATOR_HH__ */
