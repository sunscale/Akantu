/**
 * @file   segment.cc
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Thu Feb 19 2015
 * @date last modification: Thu Feb 19 2015
 *
 * @brief  CGALSegment primitive
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

#include "cgal_segment.hh"

__BEGIN_AKANTU__

template<UInt d, typename T>
CGALSegment<d, T>::CGALSegment() :
  segment(),
  source(),
  target()
{}

template<UInt d, typename T>
CGALSegment<d, T>::CGALSegment(const Point<d, T> & a, const Point<d, T> & b) :
  segment(a.getCGALPointInstance(), b.getCGALPointInstance()),
  source(a),
  target(b)
{}

template<UInt d, typename T>
T CGALSegment<d, T>::squaredLength() const {
  return segment.squared_length();
}

template<UInt d, typename T>
Segment<d, T> CGALSegment<d, T>::opposite() const {
  return CGALSegment<d, T>(segment.opposite());
}

template<UInt d, typename T>
bool CGALSegment<d, T>::operator==(const CGALSegment<d, T> & other) const {
  return this->segment == other.segment;
}

__END_AKANTU__
