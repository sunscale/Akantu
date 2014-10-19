/**
 * @file   dumper_type_traits.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Tue Sep 02 2014
 * @date last modification: Tue Sep 02 2014
 *
 * @brief  Type traits for field properties
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

#ifndef __AKANTU_DUMPER_TYPE_TRAITS_HH__
#define __AKANTU_DUMPER_TYPE_TRAITS_HH__
/* -------------------------------------------------------------------------- */
#include "element_type_map.hh"
#include "element_type_map_filter.hh"
/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__
__BEGIN_AKANTU_DUMPER__
/* -------------------------------------------------------------------------- */

template <class data, class ret, class field>
struct TypeTraits {

  //! the stored data (real, int, uint, ...)
  typedef data  data_type;
  //! the type returned by the operator *
  typedef ret   return_type;
  //! the field type (ElementTypeMap or ElementTypeMapFilter)
  typedef field field_type;
  //! the type over which we iterate
  typedef typename field_type::type it_type;
  //! the type of array (Array<T> or ArrayFilter<T>)
  typedef typename field_type::array_type array_type;
  //! the iterator over the array
  typedef typename array_type::const_vector_iterator array_iterator;
};

/* -------------------------------------------------------------------------- */

template <class T, template <class> class ret, bool filtered>
struct SingleType
  : public TypeTraits<T,ret<T>,ElementTypeMapArray<T> >{
  
};

/* -------------------------------------------------------------------------- */
template <class T, template <class> class ret>
struct SingleType<T,ret,true> : 
  public TypeTraits<T,ret<T>, ElementTypeMapArrayFilter<T> >{
  
};
/* -------------------------------------------------------------------------- */
template <class it_type, class data_type, template <class> class ret, 
	  bool filtered>
struct DualType
  : public TypeTraits<data_type,ret<data_type>, ElementTypeMapArray<it_type> >{
};

/* -------------------------------------------------------------------------- */
template <class it_type, class data_type,template <class> class ret>
struct DualType<it_type,data_type,ret,true> : 
  public TypeTraits<data_type,ret<data_type>, ElementTypeMapArrayFilter<it_type> >{
  
};
/* -------------------------------------------------------------------------- */
__END_AKANTU_DUMPER__
__END_AKANTU__

/* -------------------------------------------------------------------------- */


#endif /* __AKANTU_DUMPER_TYPE_TRAITS_HH__ */
