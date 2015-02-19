/**
 * @file   segment.hh
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

#ifndef __AKANTU_SEGMENT_HH__
#define __AKANTU_SEGMENT_HH__

#include "aka_common.hh"
#include "geometrical_primitive.hh"
#include "point.hh"

__BEGIN_AKANTU__

template<UInt d, typename T>
class Segment : public GeometricalPrimitive {

public:
  /// Default constructor
  Segment();

  /// Construct segment from two points
  Segment(const Point<d, T> & a, const Point<d, T> & b);

public:
  /// Returns first point of segment
  virtual const Point<d, T> & getSource() const;

  /// Returns last point of segment
  virtual const Point<d, T> & getTarget() const;

  /// Returns the squared length of the segment
  virtual T squaredLength() const;

  /// Returns the opposite segment
  virtual Segment<d, T> opposite() const;

public:
  /// Equality operator
  virtual bool operator==(const Segment<d, T> & other) const;

  /// Inequality operator
  virtual bool operator!=(const Segment<d, T> & other) const;
};

__END_AKANTU__

#endif // __AKANTU_SEGMENT_HH__
