/**
 * @file   triangle.cc
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Thu Feb 19 2015
 * @date last modification: Thu Feb 19 2015
 *
 * @brief  Triangle primitive
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

#include "triangle.hh"

__BEGIN_AKANTU__


template<UInt d, typename T>
Triangle<d, T>::Triangle() {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

template<UInt d, typename T>
Triangle<d, T>::Triangle(const Point<d, T> & a, const Point<d, T> & b, const Point<d, T> & c) {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

template<UInt d, typename T>
bool Triangle<d, T>::hasOnBoundary(const Point<d, T> & point) const {
  AKANTU_DEBUG_TO_IMPLEMENT();
}


template<UInt d, typename T>
bool Triangle<d, T>::hasOnBoundedSide(const Point<d, T> & point) const {
  AKANTU_DEBUG_TO_IMPLEMENT();
}


template<UInt d, typename T>
bool Triangle<d, T>::hasOnUnboundedSide(const Point<d, T> & point) const {
  AKANTU_DEBUG_TO_IMPLEMENT();
}


template<UInt d, typename T>
bool Triangle<d, T>::operator==(const Triangle<d, T> & other) const {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

template<UInt d, typename T>
bool Triangle<d, T>::operator!=(const Triangle<d, T> & other) const {
  return !(*this == other);
}

__END_AKANTU__

