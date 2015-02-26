/**
 * @file   segment.cc
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Thu Feb 19 2015
 * @date last modification: Thu Feb 19 2015
 *
 * @brief  Segment primitive
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

/* -------------------------------------------------------------------------- */

#include "segment.hh"

__BEGIN_AKANTU__

template<UInt d, typename T>
Segment<d, T>::Segment() {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

template<UInt d, typename T>
Segment<d, T>::Segment(const Point<d, T> & a, const Point<d, T> & b) {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

template<UInt d, typename T>
const Point<d, T> & Segment<d, T>::getSource() const {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

template<UInt d, typename T>
const Point<d, T> & Segment<d, T>::getTarget() const {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

template<UInt d, typename T>
T Segment<d, T>::squaredLength() const {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

template<UInt d, typename T>
Segment<d, T> Segment<d, T>::opposite() const {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

template<UInt d, typename T>
bool Segment<d, T>::operator==(const Segment<d, T> & other) const {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

template<UInt d, typename T>
bool Segment<d, T>::operator!=(const Segment<d, T> & other) const {
  return !(*this == other);
}

__END_AKANTU__
