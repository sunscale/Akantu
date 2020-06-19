/**
 * @file   aka_circular_array_inline_impl.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Fri Nov 11 2011
 * @date last modification: Wed Feb 03 2016
 *
 * @brief  implementation of circular array
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_circular_array.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class T>
inline typename CircularArray<T>::reference CircularArray<T>::
operator()(UInt i, UInt j) {
  AKANTU_DEBUG_ASSERT(end_position != start_position,
                      "The array \"" << this->id << "\" is empty");
  AKANTU_DEBUG_ASSERT(
      (i < (end_position - start_position + this->allocated_size) %
                   this->allocated_size +
               1) &&
          (j < this->nb_component),
      "The value at position [" << i << "," << j
                                << "] is out of range in array \"" << this->id
                                << "\"");
  return this->values[((i + start_position) % this->allocated_size) *
                          this->nb_component +
                      j];
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline typename CircularArray<T>::const_reference CircularArray<T>::
operator()(UInt i, UInt j) const {
  AKANTU_DEBUG_ASSERT(end_position != start_position,
                      "The array \"" << this->id << "\" is empty");
  AKANTU_DEBUG_ASSERT(
      (i < (end_position - start_position + this->allocated_size) %
                   this->allocated_size +
               1) &&
          (j < this->nb_component),
      "The value at position [" << i << "," << j
                                << "] is out of range in array \"" << this->id
                                << "\"");
  return this->values[((i + start_position) % this->allocated_size) *
                          this->nb_component +
                      j];
}

/* -------------------------------------------------------------------------- */
template <class T> inline void CircularArray<T>::makeStep() {
  AKANTU_DEBUG_IN();

  start_position = (start_position + 1) % this->allocated_size;
  end_position = (end_position + 1) % this->allocated_size;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class T>
void CircularArray<T>::printself(std::ostream & stream, int indent) const {
  std::string space;
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
    ;

  stream << space << "CircularArray<" << debug::demangle(typeid(T).name())
         << "> [" << std::endl;
  stream << space << " + start_position : " << this->start_position
         << std::endl;
  stream << space << " + end_position   : " << this->end_position << std::endl;
  Array<T>::printself(stream, indent + 1);

  stream << space << "]" << std::endl;
}

}
