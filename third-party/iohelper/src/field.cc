/**
 * @file   field.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Wed Dec 14 2011
 * @date last modification: Thu Nov 01 2012
 *
 * @brief  implementation of the Field class
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * IOHelper is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * IOHelper is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with IOHelper. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "field.hh"
/* -------------------------------------------------------------------------- */


namespace iohelper {

/* -------------------------------------------------------------------------- */
template <>
FieldType Field<Real>::getFieldDataType(){
  return REAL;
}

/* -------------------------------------------------------------------------- */
template <>
FieldType Field<UInt>::getFieldDataType(){
  return UINT;
}

/* -------------------------------------------------------------------------- */
template <>
FieldType Field<bool>::getFieldDataType(){
  return BOOL;
}

/* -------------------------------------------------------------------------- */
template <>
FieldType Field<int>::getFieldDataType(){
  return INT;
}

/* -------------------------------------------------------------------------- */
template <>
FieldType Field<long int>::getFieldDataType(){
  return INT;
}


/* -------------------------------------------------------------------------- */

}
