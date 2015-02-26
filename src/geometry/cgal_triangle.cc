/**
 * @file   cgal_triangle.cc
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Thu Feb 19 2015
 * @date last modification: Thu Feb 19 2015
 *
 * @brief  CGALTriangle primitive using CGAL implementation
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

#include "cgal_triangle.hh"

__BEGIN_AKANTU__


template<UInt d, typename T>
CGALTriangle<d, T>::CGALTriangle() :
  triangle()
{}

template<UInt d, typename T>
CGALTriangle<d, T>::CGALTriangle(const Point<d, T> & a, const Point<d, T> & b, const Point<d, T> & c) :
  triangle(a.getCGALPointInstance(), b.getCGALPointInstance(), c.getCGALPointInstance())
{}

template<UInt d, typename T>
bool CGALTriangle<d, T>::hasOnBoundary(const Point<d, T> & point) const {
  return triangle.has_on_boundary(point.getCGALPointInstance());
}


template<UInt d, typename T>
bool CGALTriangle<d, T>::hasOnBoundedSide(const Point<d, T> & point) const {
  return triangle.has_on_bounded_side(point.getCGALPointInstance());
}


template<UInt d, typename T>
bool CGALTriangle<d, T>::hasOnUnboundedSide(const Point<d, T> & point) const {
  return triangle.has_on_unbounded_side(point.getCGALPointInstance());
}


template<UInt d, typename T>
bool CGALTriangle<d, T>::operator==(const CGALTriangle<d, T> & other) const {
  return this->triangle == other.triangle;
}

__END_AKANTU__

