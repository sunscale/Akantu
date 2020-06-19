/**
 * @file   aka_array.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Implementation of akantu::Array
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
#include <memory>
#include <utility>

/* -------------------------------------------------------------------------- */
#include "aka_array.hh"
#include "aka_common.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
/* Functions ArrayBase                                                       */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template <> UInt Array<Real>::find(const Real & elem) const {
  AKANTU_DEBUG_IN();

  Real epsilon = std::numeric_limits<Real>::epsilon();
  auto it = std::find_if(begin(), end(), [&elem, &epsilon](auto && a) {
    return std::abs(a - elem) <= epsilon;
  });

  AKANTU_DEBUG_OUT();
  return (it != end()) ? end() - it : UInt(-1);
}

/* -------------------------------------------------------------------------- */
template <>
Array<ElementType> & Array<ElementType>::operator*=(__attribute__((unused))
                                                    const ElementType & alpha) {
  AKANTU_TO_IMPLEMENT();
  return *this;
}

template <>
Array<ElementType> & Array<ElementType>::
operator-=(__attribute__((unused)) const Array<ElementType> & vect) {
  AKANTU_TO_IMPLEMENT();
  return *this;
}

template <>
Array<ElementType> & Array<ElementType>::
operator+=(__attribute__((unused)) const Array<ElementType> & vect) {
  AKANTU_TO_IMPLEMENT();
  return *this;
}

template <>
Array<char> & Array<char>::operator*=(__attribute__((unused))
                                      const char & alpha) {
  AKANTU_TO_IMPLEMENT();
  return *this;
}

template <>
Array<char> & Array<char>::operator-=(__attribute__((unused))
                                      const Array<char> & vect) {
  AKANTU_TO_IMPLEMENT();
  return *this;
}

template <>
Array<char> & Array<char>::operator+=(__attribute__((unused))
                                      const Array<char> & vect) {
  AKANTU_TO_IMPLEMENT();
  return *this;
}

} // namespace akantu
