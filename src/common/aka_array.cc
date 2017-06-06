/**
 * @file   aka_array.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Aug 18 2015
 *
 * @brief  Implementation of akantu::Array
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include <memory>
#include <utility>

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_array.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/* Functions ArrayBase                                                       */
/* -------------------------------------------------------------------------- */
ArrayBase::ArrayBase(ID id)
    : id(std::move(id)), allocated_size(0), size(0), nb_component(1),
      size_of_type(0) {}

/* -------------------------------------------------------------------------- */
ArrayBase::~ArrayBase() = default;

/* -------------------------------------------------------------------------- */
void ArrayBase::printself(std::ostream & stream, int indent) const {
  std::string space;
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
    ;
  stream << space << "ArrayBase [" << std::endl;
  stream << space << " + size             : " << size << std::endl;
  stream << space << " + nb component     : " << nb_component << std::endl;
  stream << space << " + allocated size   : " << allocated_size << std::endl;
  Real mem_size = (allocated_size * nb_component * size_of_type) / 1024.;
  stream << space << " + size of type     : " << size_of_type << "B"
         << std::endl;
  stream << space << " + memory allocated : " << mem_size << "kB" << std::endl;
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
template <> Int Array<Real>::find(const Real & elem) const {
  AKANTU_DEBUG_IN();
  UInt i = 0;
  Real epsilon = std::numeric_limits<Real>::epsilon();
  for (; (i < size) && (fabs(values[i] - elem) <= epsilon); ++i)
    ;

  AKANTU_DEBUG_OUT();
  return (i == size) ? -1 : (Int)i;
}

/* -------------------------------------------------------------------------- */
template <>
Array<ElementType> & Array<ElementType>::operator*=(__attribute__((unused))
                                                    const ElementType & alpha) {
  AKANTU_DEBUG_TO_IMPLEMENT();
  return *this;
}

template <>
Array<ElementType> & Array<ElementType>::
operator-=(__attribute__((unused)) const Array<ElementType> & vect) {
  AKANTU_DEBUG_TO_IMPLEMENT();
  return *this;
}

template <>
Array<ElementType> & Array<ElementType>::
operator+=(__attribute__((unused)) const Array<ElementType> & vect) {
  AKANTU_DEBUG_TO_IMPLEMENT();
  return *this;
}

template <>
Array<char> & Array<char>::operator*=(__attribute__((unused))
                                      const char & alpha) {
  AKANTU_DEBUG_TO_IMPLEMENT();
  return *this;
}

template <>
Array<char> & Array<char>::operator-=(__attribute__((unused))
                                      const Array<char> & vect) {
  AKANTU_DEBUG_TO_IMPLEMENT();
  return *this;
}

template <>
Array<char> & Array<char>::operator+=(__attribute__((unused))
                                      const Array<char> & vect) {
  AKANTU_DEBUG_TO_IMPLEMENT();
  return *this;
}

__END_AKANTU__
