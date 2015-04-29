/**
 * @file   geom_helper_functions.hh
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Wed Mar 4 2015
 * @date last modification: Thu Mar 5 2015
 *
 * @brief  Helper functions for the computational geometry algorithms
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

#ifndef _AKANTU_GEOM_HELPER_FUNCTIONS_HH__
#define _AKANTU_GEOM_HELPER_FUNCTIONS_HH__

#include "aka_common.hh"
#include "aka_math.hh"
#include "tree_type_helper.hh"

#include <CGAL/Cartesian.h>

__BEGIN_AKANTU__

typedef CGAL::Cartesian<Real> K;

#define EPS 1e-10

/// Fuzzy compare of two points
inline bool comparePoints(const K::Point_3 & a, const K::Point_3 & b) {
  Math::setTolerance(EPS);
  return Math::are_float_equal(a.x(), b.x()) && Math::are_float_equal(a.y(), b.y());
}

/// Fuzzy compare of two segments
inline bool compareSegments(const K::Segment_3 & a, const K::Segment_3 & b) {
  return (comparePoints(a.source(), b.source()) && comparePoints(a.target(), b.target())) ||
         (comparePoints(a.source(), b.target()) && comparePoints(a.target(), b.source()));
}

/// Compare segment pairs
inline bool compareSegmentPairs(const std::pair<K::Segment_3, UInt> & a, const std::pair<K::Segment_3, UInt> & b) {
  return compareSegments(a.first, b.first);
}

/// Pair ordering operator based on second member
inline bool comparePairElement(const std::pair<K::Segment_3, UInt> & a, const std::pair<K::Segment_3, UInt> & b) {
  return a.second < b.second;
}

/// Pair ordering operator based on first member
inline bool lessSegmentPair(const std::pair<K::Segment_3, UInt> & a, const std::pair<K::Segment_3, UInt> & b) {
  return CGAL::compare_lexicographically(a.first.min(), b.first.min()) || CGAL::compare_lexicographically(a.first.max(), b.first.max());
}

/* -------------------------------------------------------------------------- */
/* Predicates                                                                 */
/* -------------------------------------------------------------------------- */

/// Predicate used to eliminate faces of mesh not belonging to a specific element
template <UInt dim, ElementType el_type>
class BelongsNotToElement {

public:
  BelongsNotToElement(UInt el):
    el(el)
  {}

  bool operator()(const typename TreeTypeHelper<dim, el_type>::primitive_type & primitive) {
    return primitive.id() != el;
  }

protected:
  const UInt el;
};

/// Predicate used to determine if point is on edge of faces of mesh
template <UInt dim, ElementType el_type>
class HasOnEdge {

public:
  HasOnEdge(const K::Point_3 & point):
    point(point)
  {}

  bool operator()(const typename TreeTypeHelper<dim, el_type>::primitive_type & primitive) {
    return primitive.has_on(point);
  }

protected:
  const K::Point_3 & point;
};

/// Predicate used to determine if two segments are equal
class IsSameSegment {

public:
  IsSameSegment(const K::Segment_3 & segment):
    segment(segment)
  {}

  bool operator()(const std::pair<K::Segment_3, UInt> & test_pair) {
    return compareSegments(segment, test_pair.first);
  }
protected:
  const K::Segment_3 segment;
};

__END_AKANTU__

#endif // _AKANTU_GEOM_HELPER_FUNCTIONS_HH__

