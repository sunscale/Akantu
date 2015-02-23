/**
 * @file   point.hh
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Thu Feb 19 2015
 * @date last modification: Thu Feb 19 2015
 *
 * @brief  Point primitive
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

#ifndef __AKANTU_POINT_HH__
#define __AKANTU_POINT_HH__

#include "aka_common.hh"
#include "geometrical_primitive.hh"

__BEGIN_AKANTU__

template<UInt d, typename T = Real>
class Point : public GeometricalPrimitive {

public:
  /// Default constructor
  Point();

  /// Construct from akantu vector
  explicit Point(const Vector<T> & vec);

  /// Copy constructor
  explicit Point(const Point<d, T> & other);

public:
  /// Equality operator
  virtual bool operator==(const Point<d, T> & other) const;

  /// Inequality operator
  virtual bool operator!=(const Point<d, T> & other) const;

  /// Returns the i-th componenent
  virtual T operator[](UInt i) const;
};

__END_AKANTU__

#endif // __AKANTU_POINT_HH__
