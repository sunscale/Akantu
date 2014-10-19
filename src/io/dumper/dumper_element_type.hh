/**
 * @file   dumper_element_type.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Tue Sep 02 2014
 * @date last modification: Tue Sep 02 2014
 *
 * @brief  ElementType Field
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

__BEGIN_AKANTU__
__BEGIN_AKANTU_DUMPER__
/* -------------------------------------------------------------------------- */

template <class types>
class element_type_field_iterator 
  : public element_iterator<types, element_type_field_iterator> {
public:
  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */ 
  typedef element_iterator<types, dumper::element_type_field_iterator> parent;

  typedef typename types::data_type   data_type;
  typedef typename types::return_type return_type;
  typedef typename types::array_iterator array_iterator;
  typedef typename types::field_type field_type;

public:

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

  element_type_field_iterator(const field_type & field,
                              const typename field_type::type_iterator & t_it,
                              const typename field_type::type_iterator & t_it_end,
                              const array_iterator & array_it,
                              const array_iterator & array_it_end,
                              const GhostType ghost_type = _not_ghost) :
    parent(field, t_it, t_it_end, array_it, array_it_end, ghost_type) { }

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
  
  return_type operator*() {
    data_type type = getIOHelperType(*this->tit);
    return return_type(1, type);
  }
};

/* -------------------------------------------------------------------------- */


template <bool filtered = false>
class ElementTypeField 
  : public GenericElementalField<DualType<UInt,iohelper::ElemType,Vector,filtered>,
 				 element_type_field_iterator> {
public:
  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */  

  typedef DualType<UInt,iohelper::ElemType,Vector,filtered> types;
  typedef element_type_field_iterator<types> iterator;
  typedef typename types::field_type field_type;
  typedef GenericElementalField<types,element_type_field_iterator> parent;

public:

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

  ElementTypeField(const field_type & field,
		   UInt spatial_dimension = _all_dimensions,
		   GhostType ghost_type = _not_ghost,
		   ElementKind element_kind = _ek_not_defined) :
    parent(field, spatial_dimension, ghost_type, element_kind) {
    this->homogeneous = true;
  }

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

  UInt getDim() { return 1; }
};


/* -------------------------------------------------------------------------- */

__END_AKANTU_DUMPER__
__END_AKANTU__
/* -------------------------------------------------------------------------- */
