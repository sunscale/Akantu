/**
 * @file   dumper_filtered_connectivity.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Tue Sep 02 2014
 * @date last modification: Tue Sep 02 2014
 *
 * @brief  FilteredConnectivities field
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

// #include "dumper_generic_elemental_field.hh"
// /* -------------------------------------------------------------------------- */
// __BEGIN_AKANTU__
// __BEGIN_AKANTU_DUMPER__
// /* -------------------------------------------------------------------------- */

// template <class types>
// class filtered_connectivity_field_iterator 
//   : public element_iterator<types, filtered_connectivity_field_iterator> {

// public:

//   /* ------------------------------------------------------------------------ */
//   /* Typedefs                                                                 */
//   /* ------------------------------------------------------------------------ */  
  

//   typedef element_iterator<types, dumper::filtered_connectivity_field_iterator> parent;
//   typedef typename types::return_type return_type;
//   typedef typename types::field_type field_type;
//   typedef typename types::array_iterator array_iterator;

// public:

//   /* ------------------------------------------------------------------------ */
//   /* Constructors/Destructors                                                 */
//   /* ------------------------------------------------------------------------ */

//   filtered_connectivity_field_iterator(const field_type & field,
// 				       const typename field_type::type_iterator & t_it,
//                                        const typename field_type::type_iterator & t_it_end,
//                                        const array_iterator & array_it,
//                                        const array_iterator & array_it_end,
// 				       const GhostType ghost_type = _not_ghost) :
//     parent(field, t_it, t_it_end, array_it, array_it_end, ghost_type) { }

//   /* ------------------------------------------------------------------------ */
//   /* Methods                                                                  */
//   /* ------------------------------------------------------------------------ */
//   bool operator!=(const iterator & it) const {
//     return (ghost_type != it.ghost_type)
//       || (tit != it.tit || ((array_it != it.array_it) || sub_element != it.sub_element) );
//   }

//   iterator & operator++() {
//     if (!this->sub_element)
//       this->sub_element += 1;
//     else {
//       ++array_it;
//       this->sub_element = 0;
//       while(array_it == array_it_end && tit != tit_end) {
// 	++tit;
// 	if(tit != tit_end) {

// 	  const array_type & vect = field(*tit, ghost_type);
// 	  UInt _nb_data_per_elem = getNbDataPerElem(*tit);
// 	  UInt nb_component = vect.getNbComponent();
// 	  UInt size = (vect.getSize() * nb_component) / _nb_data_per_elem;

// 	  array_it       = vect.begin_reinterpret(_nb_data_per_elem,size);
// 	  array_it_end   = vect.end_reinterpret  (_nb_data_per_elem,size);
// 	}
//       }
//     }
//   }

//   return_type operator*(){
//     const Vector<UInt> & element_connect = *this->array_it;
//     switch (sub_element) {
//     case 0:
//       UInt * sub_connec_ptr = InterpolationElement<ElementClassProperty<*tit>::interpolation_type>::sub_element_connectivity[sub_element];
//       UInt nb_nodes_sub_el = ElementClass<ElementClassProperty<*tit>::sub_element_type_1>::getNbNodesPerInterpolationElement(); break;
//     case 1:
//       UInt * sub_connec_ptr = InterpolationElement<ElementClassProperty<*tit>::interpolation_type>::sub_element_connectivity[sub_element];
//       UInt nb_nodes_sub_el = ElementClass<ElementClassProperty<*tit>::sub_element_type_2>::getNbNodesPerInterpolationElement(); break;
//     }

//     const Vector<UInt> sub_element_connect(sub_connect_ptr, nb_nodes_sub_el);
//     return sub_element_connect;
//   }


//   /* ------------------------------------------------------------------------ */
//   /* Class Members                                                            */
//   /* ------------------------------------------------------------------------ */
  
// private:
//   UInt sub_element;
// };

// /* -------------------------------------------------------------------------- */

// class FilteredConnectivityField : 
//   public GenericElementalField<SingleType<UInt,Vector,true>,
//  			       filtered_connectivity_field_iterator> {

//   /* ------------------------------------------------------------------------ */
//   /* Typedefs                                                                 */
//   /* ------------------------------------------------------------------------ */  

// public:

//   typedef SingleType<UInt,Vector,true> types;
//   typedef filtered_connectivity_field_iterator<types> iterator;
//   typedef types::field_type field_type;
//   typedef GenericElementalField<types,filtered_connectivity_field_iterator> parent;

//   /* ------------------------------------------------------------------------ */
//   /* Constructors/Destructors                                                 */
//   /* ------------------------------------------------------------------------ */

//   FilteredConnectivityField(const field_type & field,
//  			    const Array<UInt> & nodal_filter,
// 			    UInt spatial_dimension = _all_dimensions,
// 			    GhostType ghost_type = _not_ghost,
// 			    ElementKind element_kind = _ek_not_defined) :
//     parent(field, spatial_dimension, ghost_type, element_kind),
//     nodal_filter(nodal_filter) { }

//   /* ------------------------------------------------------------------------ */
//   /* Methods                                                                  */
//   /* ------------------------------------------------------------------------ */
  
//   iterator begin() {
//     iterator it = parent::begin();
//     it.setNodalFilter(nodal_filter);
//     return it;
//   }

//   iterator end() {
//     iterator it = parent::end();
//     it.setNodalFilter(nodal_filter);
//     return it;
//   }

//   /* ------------------------------------------------------------------------ */
//   /* Class Members                                                            */
//   /* ------------------------------------------------------------------------ */


// private:
//   const Array<UInt> & nodal_filter;
// };


// /* -------------------------------------------------------------------------- */

// __END_AKANTU_DUMPER__
// __END_AKANTU__

// /* -------------------------------------------------------------------------- */
