/**
 * @file   cgal_point.hh
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

#ifndef __AKANTU_CGAL_POINT_HH__
#define __AKANTU_CGAL_POINT_HH__

#include "aka_common.hh"
#include "point.hh"

#include <CGAL/Cartesian.h>

#define COMMA ,

__BEGIN_AKANTU__

/// Definition of the type for the CGAL point instance in class CGALPoint
template<UInt d, typename T> class CGAL_point_instance;

/// Special case for 1D
template<typename T> class CGAL_point_instance<1, T> {
  
public:
  CGAL_point_instance() { AKANTU_DEBUG_TO_IMPLEMENT(); }
};

/// 2D and 3D use CGAL objects
template<typename T> class CGAL_point_instance<2, T> : public CGAL::Cartesian<T>::Point_2 {};
template<typename T> class CGAL_point_instance<3, T> : public CGAL::Cartesian<T>::Point_3 {};

/* -------------------------------------------------------------------------- */

template<UInt d, typename T = Real>
class CGALPoint : public Point<d, T> {

public:
  /// Default constructor, constructs a point on origin
  CGALPoint();

  /// Constructs a point from vector
  explicit CGALPoint(const Vector<T> & vector);

  /// Copy constructor
  explicit CGALPoint(const Point<d, T> & other);

public:
  /// Getter for point
  AKANTU_GET_MACRO(CGALPointInstance, point, CGAL_point_instance<d COMMA T>);

public:
  /// Equality operator
  bool operator==(const CGALPoint<d, T> & other) const;

  /// Inequality operator
  bool operator!=(const CGALPoint<d, T> & other) const;

  /// Returns the i-th component
  T operator[](UInt i) const;

protected:
  CGAL_point_instance<d, T> point;
};

__END_AKANTU__

#endif // __AKANTU_CGAL_POINT_HH__
