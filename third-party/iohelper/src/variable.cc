/**
 * @file   variable.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Tue Jun 04 2013
 * @date last modification: Tue Jun 04 2013
 *
 * @brief  implementation of the Variable class
 *
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "variable.hh"
/* -------------------------------------------------------------------------- */

namespace iohelper {

/* -------------------------------------------------------------------------- */
template <>
VarType Variable<Real>::getVarDataType(){
  return REAL;
}

/* -------------------------------------------------------------------------- */
template <>
VarType Variable<UInt>::getVarDataType(){
  return UINT;
}

/* -------------------------------------------------------------------------- */
template <>
VarType Variable<bool>::getVarDataType(){
  return BOOL;
}

/* -------------------------------------------------------------------------- */
template <>
VarType Variable<Int>::getVarDataType(){
  return INT;
}

/* -------------------------------------------------------------------------- */

}
