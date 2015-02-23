/**
 * @file   cgal_segment.hh
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Thu Feb 19 2015
 * @date last modification: Thu Feb 19 2015
 *
 * @brief  Segment primitive using CGAL implementation
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

#ifndef __AKANTU_CGAL_SEGMENT_HH__
#define __AKANTU_CGAL_SEGMENT_HH__

#include "aka_common.hh"
#include "segment.hh"

#include <CGAL/Cartesian.h>

// TODO : find less ugly solution for the macro-template problem
#define COMMA ,

__BEGIN_AKANTU__

/// Definition of the type for the CGAL segment instance in CGALSegment
template<UInt d, typename T> class CGAL_segment_instance;

/// Special case for 1D
template<typename T> class CGAL_segment_instance<1, T> {
public:
  CGAL_segment_instance() {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }
};

/// 2D and 3D use CGAL objects
template<typename T> class CGAL_segment_instance<2, T> : public CGAL::Cartesian<T>::Segment_2 {};
template<typename T> class CGAL_segment_instance<3, T> : public CGAL::Cartesian<T>::Segment_3 {};

/* -------------------------------------------------------------------------- */


template<UInt d, typename T = Real>
class CGALSegment : public Segment<d, T> {

public:
  /// Default constructor
  CGALSegment();

  /// Construct segment from two points
  CGALSegment(const Point<d, T> & a, const Point<d, T> & b);

public:
  /// Returns the squared length of the segment
  virtual T squaredLength() const;

  /// Returns the opposite segment
  virtual Segment<d, T> opposite() const;

public:
  /// Getter for CGAL instance
  AKANTU_GET_MACRO(CGALSegmentInstance, segment, CGAL_segment_instance<d COMMA T>);

  /// Getter for source
  AKANTU_GET_MACRO(Source, source, Point<d COMMA T>);

  /// Getter for target
  AKANTU_GET_MACRO(Target, target, Point<d COMMA T>);

public:
  /// Equality operator
  virtual bool operator==(const CGALSegment<d, T> & other) const;

protected:
  CGAL_segment_instance<d, T> segment;

  Point<d, T> source, target;
};

__END_AKANTU__

#endif // __AKANTU_CGAL_SEGMENT_HH__
