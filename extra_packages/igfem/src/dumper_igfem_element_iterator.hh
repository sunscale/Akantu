/**
 * @file   dumper_igfem_element_iterator.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  Iterators for IGFEM elemental fields
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

#ifndef AKANTU_DUMPER_IGFEM_ELEMENT_ITERATOR_HH_
#define AKANTU_DUMPER_IGFEM_ELEMENT_ITERATOR_HH_
/* -------------------------------------------------------------------------- */
#include "element.hh"
#include "igfem_helper.hh"
/* -------------------------------------------------------------------------- */
namespace akantu {
namespace dumpers {
/* -------------------------------------------------------------------------- */

template <class types, template <class> class final_iterator>
class igfem_element_iterator {
  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  typedef typename types::it_type it_type;
  typedef typename types::field_type field_type;
  typedef typename types::array_type array_type;
  typedef typename types::array_iterator array_iterator;
  typedef final_iterator<types> iterator;

public:
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

  igfem_element_iterator(const field_type & field,
                         const typename field_type::type_iterator & t_it,
                         const typename field_type::type_iterator & t_it_end,
                         const array_iterator & array_it,
                         const array_iterator & array_it_end,
                         const GhostType ghost_type = _not_ghost,
                         UInt sub_element = 0)
      : field(field), tit(t_it), tit_end(t_it_end), array_it(array_it),
        array_it_end(array_it_end), ghost_type(ghost_type),
        sub_element(sub_element) {}

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

public:
  bool operator!=(const iterator & it) const {
    return (ghost_type != it.ghost_type) ||
           (tit != it.tit ||
            ((array_it != it.array_it) || sub_element != it.sub_element));
  }

  iterator & operator++() {
    if (!this->sub_element)
      this->sub_element += 1;
    else {
      ++array_it;
      this->sub_element = 0;
      while (array_it == array_it_end && tit != tit_end) {
        ++tit;
        if (tit != tit_end) {
          const array_type & vect = field(*tit, ghost_type);
          UInt _nb_data_per_elem = getNbDataPerElem(*tit);
          UInt nb_component = vect.getNbComponent();
          UInt size = (vect.getSize() * nb_component) / _nb_data_per_elem;

          array_it = vect.begin_reinterpret(_nb_data_per_elem, size);
          array_it_end = vect.end_reinterpret(_nb_data_per_elem, size);
        }
      }
    }
    return *(static_cast<iterator *>(this));
  }

  ElementType getType() {
    ElementType sub_type = IGFEMHelper::getSubElementType(*tit, sub_element);
    return sub_type;
  }

  /// get IOHelperType for sub-element
  UInt element_type() { return getIOHelperType(this->getType()); }

  /// get current parent element????
  Element getCurrentElement() {
    return Element(*tit, array_it.getCurrentIndex());
  }

  UInt getNbDataPerElem(ElementType type) const {
    /// nb of data per parent element!
    if (!nb_data_per_elem.exists(type, ghost_type))
      return field(type, ghost_type).getNbComponent();

    return nb_data_per_elem(type, ghost_type);
  }

  void setNbDataPerElem(const ElementTypeMap<UInt> & nb_data) {
    /// nb of data per parent element!
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
  /// index of sub-element
  UInt sub_element;
  /// sub_element end
  UInt sub_element_end;
};

/* -------------------------------------------------------------------------- */
template <typename types>
class igfem_elemental_field_iterator
    : public igfem_element_iterator<types, igfem_elemental_field_iterator> {
public:
  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */

  typedef igfem_element_iterator<
      types, ::akantu::dumpers::igfem_elemental_field_iterator>
      parent;
  typedef typename types::it_type it_type;
  typedef typename types::return_type return_type;
  typedef typename types::field_type field_type;
  typedef typename types::array_iterator array_iterator;

public:
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

  igfem_elemental_field_iterator(
      const field_type & field, const typename field_type::type_iterator & t_it,
      const typename field_type::type_iterator & t_it_end,
      const array_iterator & array_it, const array_iterator & array_it_end,
      const GhostType ghost_type = _not_ghost, UInt sub_element = 0)
      : parent(field, t_it, t_it_end, array_it, array_it_end, ghost_type,
               sub_element) {}

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

  return_type operator*() { return *this->array_it; }

private:
};

/* -------------------------------------------------------------------------- */
} // namespace dumpers
} // namespace akantu
/* -------------------------------------------------------------------------- */

#endif /* AKANTU_DUMPER_IGFEM_ELEMENT_ITERATOR_HH_ */
