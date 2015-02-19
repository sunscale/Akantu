/**
 * @file   cgal_point.cc
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Thu Feb 19 2015
 * @date last modification: Thu Feb 19 2015
 *
 * @brief  Point primitive using CGAL implementation
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

#include "cgal_point.hh"

__BEGIN_AKANTU__

template<UInt d, typename T>
CGALPoint<d, T>::CGALPoint() :
  point()
{}

template<UInt d, typename T>
CGALPoint<d, T>::CGALPoint(const Vector<T> & vector) :
  point()
{
  for (UInt i = 0 ; i < d ; i++) {
    point[i] = vector[i];
  }
}

template<UInt d, typename T>
CGALPoint<d, T>::CGALPoint(const Point<d, T> & other) :
  point()
{
  for (UInt i = 0 ; i < d ; i++) {
    point[i] = other[i];
  }
}

template<UInt d, typename T>
bool CGALPoint<d, T>::operator==(const CGALPoint<d, T> & other) const {
  return this->point == other.point;
}

template<UInt d, typename T>
T CGALPoint<d, T>::operator[](UInt i) const {
  return point[i];
}

__END_AKANTU__
