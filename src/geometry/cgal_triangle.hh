/**
 * @file   cgal_triangle.hh
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Thu Feb 19 2015
 * @date last modification: Thu Feb 19 2015
 *
 * @brief  Triangle primitive using CGAL implementation
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

#ifndef __AKANTU_CGAL_TRIANGLE_HH__
#define __AKANTU_CGAL_TRIANGLE_HH__

#include "triangle.hh"
#include "point.hh"

#include <CGAL/Cartesian.h>

__BEGIN_AKANTU__

/// Definition of the type for the CGAL triangle instance in CGALTriangle
template<UInt d, typename T> class CGAL_triangle_instance;

/// Special case for 1D
template<typename T> class CGAL_triangle_instance<1, T> {
public:
  CGAL_triangle_instance() {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }
};

/// 2D and 3D use CGAL objects
template<typename T> class CGAL_triangle_instance<2, T> : public CGAL::Cartesian<T>::Triangle_2 {};
template<typename T> class CGAL_triangle_instance<3, T> : public CGAL::Cartesian<T>::Triangle_3 {};

/* -------------------------------------------------------------------------- */

template<UInt d, typename T>
class CGALTriangle : public Triangle<d, T> {

public:
  /// Default constructor
  CGALTriangle();

  /// Construct triangle from 3 points
  CGALTriangle(const Point<d, T> & a, const Point<d, T> & b, const Point<d, T> & c);

public:
  /// Check if a point is on the triangle boundary
  virtual bool hasOnBoundary(const Point<d, T> & point) const;

  /// Check if the point is inside the triangle
  virtual bool hasOnBoundedSide(const Point<d, T> & point) const;

  /// Check if the point is outside the triangle
  virtual bool hasOnUnboundedSide(const Point<d, T> & point) const;

public:
  /// Equality operator
  virtual bool operator==(const CGALTriangle<d, T> & other) const;

protected:
  CGAL_triangle_instance<d, T> triangle;

};

__END_AKANTU__

#endif // __AKANTU_CGAL_TRIANGLE_HH__
